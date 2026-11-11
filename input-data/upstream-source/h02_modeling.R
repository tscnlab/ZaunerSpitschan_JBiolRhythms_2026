# H02 temporal-model fitting, diagnostics, predictions, and variation summaries.

h02_prepare_fit_data <- function(frame) {
  site_levels <- sort(unique(frame$site))
  participant_levels <- sort(unique(frame$participant_key))
  day_levels <- sort(unique(frame$participant_day_key))
  frame |>
    dplyr::mutate(
      site = factor(.data$site, levels = site_levels),
      site_smooth = ordered(
        as.character(.data$site),
        levels = site_levels
      ),
      participant = factor(
        .data$participant_key,
        levels = participant_levels
      ),
      participant_day = factor(
        .data$participant_day_key,
        levels = day_levels
      )
    ) |>
    dplyr::arrange(
      .data$site,
      .data$participant,
      .data$local_date,
      .data$source_utc_start,
      .data$clock_bin
    )
}

h02_fit_bam <- function(
  formula,
  data,
  method,
  rho = 0,
  discrete = TRUE
) {
  warnings <- character()
  environment(formula) <- environment()
  fit <- withCallingHandlers(
    mgcv::bam(
      formula = formula,
      data = data,
      method = method,
      discrete = discrete,
      rho = rho,
      AR.start = AR_start,
      knots = list(time_hour = c(0, 24)),
      nthreads = 1L,
      gc.level = 1L,
      drop.unused.levels = TRUE
    ),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  attr(fit, "h02_warnings") <- unique(warnings)
  fit
}

h02_boundary_lag_correlation <- function(residual, AR_start, lag = 1L) {
  if (lag < 1L || lag != as.integer(lag)) {
    h02_abort("lag must be a positive integer")
  }
  sequence_id <- cumsum(AR_start)
  index <- seq_along(residual)
  earlier <- index - lag
  eligible <- earlier >= 1L
  eligible[eligible] <- sequence_id[index[eligible]] ==
    sequence_id[earlier[eligible]]
  x <- residual[index[eligible]]
  y <- residual[earlier[eligible]]
  complete <- is.finite(x) & is.finite(y)
  if (sum(complete) < 3L) {
    return(c(correlation = NA_real_, pairs = sum(complete)))
  }
  c(
    correlation = stats::cor(x[complete], y[complete]),
    pairs = sum(complete)
  )
}

h02_estimate_rho <- function(preliminary_fit, data) {
  raw <- stats::residuals(preliminary_fit, type = "response")
  estimate <- unname(
    h02_boundary_lag_correlation(raw, data$AR_start, lag = 1L)[
      "correlation"
    ]
  )
  if (!is.finite(estimate)) {
    h02_abort("Could not estimate a boundary-aware lag-1 correlation")
  }
  pmin(0.95, pmax(-0.95, estimate))
}

h02_model_row <- function(fit, model_id, run_id, rho) {
  ll <- stats::logLik(fit)
  convergence <- if (
    is.list(fit$outer.info) &&
      !is.null(fit$outer.info$conv)
  ) {
    as.character(fit$outer.info$conv)
  } else if (
    is.list(fit$mgcv.conv) &&
      !is.null(fit$mgcv.conv$fully.converged)
  ) {
    if (isTRUE(fit$mgcv.conv$fully.converged)) {
      "full convergence"
    } else {
      "not fully converged"
    }
  } else if (!is.null(fit$mgcv.conv)) {
    paste(fit$mgcv.conv, collapse = "; ")
  } else {
    "not reported"
  }
  tibble::tibble(
    run_id = run_id,
    model_id = model_id,
    method = as.character(fit$method),
    n = stats::nobs(fit),
    log_likelihood = as.numeric(ll),
    likelihood_df = as.numeric(attr(ll, "df")),
    AIC = stats::AIC(fit),
    total_edf = sum(fit$edf),
    rank = fit$rank,
    coefficients = length(stats::coef(fit)),
    scale = fit$sig2,
    rho = rho,
    convergence = convergence,
    warnings = paste(attr(fit, "h02_warnings"), collapse = " | ")
  )
}

h02_compare_models <- function(
  reduced,
  full,
  reduced_id,
  full_id,
  comparison_id,
  run_id
) {
  ll_reduced <- stats::logLik(reduced)
  ll_full <- stats::logLik(full)
  statistic <- 2 * (as.numeric(ll_full) - as.numeric(ll_reduced))
  df <- as.numeric(attr(ll_full, "df") - attr(ll_reduced, "df"))
  log_p <- if (is.finite(statistic) && statistic >= 0 && df > 0) {
    stats::pchisq(
      statistic,
      df = df,
      lower.tail = FALSE,
      log.p = TRUE
    )
  } else {
    NA_real_
  }
  p <- exp(log_p)
  tibble::tibble(
    run_id = run_id,
    comparison_id = comparison_id,
    reduced_model = reduced_id,
    full_model = full_id,
    reduced_AIC = stats::AIC(reduced),
    full_AIC = stats::AIC(full),
    delta_AIC_reduced_minus_full = stats::AIC(reduced) - stats::AIC(full),
    test_method = "restricted-likelihood diagnostic from fitted log likelihoods",
    test_statistic = statistic,
    test_df = df,
    p_raw = p,
    log_p_raw = log_p
  )
}

h02_joint_smooth_wald <- function(fit, label_pattern) {
  smooth_indices <- which(grepl(
    label_pattern,
    h02_smooth_labels(fit),
    fixed = TRUE
  ))
  if (length(smooth_indices) == 0L) {
    h02_abort("No smooth labels matched '%s'", label_pattern)
  }
  coefficient_indices <- unique(unlist(lapply(
    fit$smooth[smooth_indices],
    function(smooth) seq.int(smooth$first.para, smooth$last.para)
  )))
  estimate <- stats::coef(fit)[coefficient_indices]
  covariance <- fit$Vp[
    coefficient_indices,
    coefficient_indices,
    drop = FALSE
  ]
  eigen_decomposition <- eigen(
    covariance,
    symmetric = TRUE,
    only.values = FALSE
  )
  tolerance <- max(eigen_decomposition$values) *
    sqrt(.Machine$double.eps)
  retained <- eigen_decomposition$values > tolerance
  rank <- sum(retained)
  if (rank < 1L) {
    h02_abort("Joint smooth covariance has zero numerical rank")
  }
  projected <- drop(
    crossprod(eigen_decomposition$vectors[, retained, drop = FALSE], estimate)
  )
  statistic <- sum(
    projected^2 / eigen_decomposition$values[retained]
  )
  log_p <- stats::pchisq(
    statistic,
    df = rank,
    lower.tail = FALSE,
    log.p = TRUE
  )
  c(
    statistic = statistic,
    df = rank,
    p = exp(log_p),
    log_p = log_p,
    coefficients = length(coefficient_indices)
  )
}

h02_fit_primary_structure <- function(frame, run_id) {
  data <- h02_prepare_fit_data(frame)
  formulas <- h02_formula_set()
  preliminary <- h02_fit_bam(
    formulas$site_pattern,
    data,
    method = "fREML",
    rho = 0
  )
  rho <- h02_estimate_rho(preliminary, data)
  candidate_ids <- c(
    "no_site",
    "site_pattern",
    "no_participant_pattern",
    "no_participant_day"
  )
  candidates <- lapply(candidate_ids, function(id) {
    message("  fitting comparable fREML structure: ", id)
    h02_fit_bam(
      formulas[[id]],
      data,
      method = "fREML",
      rho = rho
    )
  })
  names(candidates) <- candidate_ids
  model_table <- dplyr::bind_rows(lapply(candidate_ids, function(id) {
    h02_model_row(candidates[[id]], id, run_id, rho)
  }))
  comparisons <- dplyr::bind_rows(
    h02_compare_models(
      candidates$no_site,
      candidates$site_pattern,
      "no_site",
      "site_pattern",
      "site_pattern_vs_no_site",
      run_id
    ),
    h02_compare_models(
      candidates$no_participant_pattern,
      candidates$site_pattern,
      "no_participant_pattern",
      "site_pattern",
      "participant_pattern",
      run_id
    ),
    h02_compare_models(
      candidates$no_participant_day,
      candidates$site_pattern,
      "no_participant_day",
      "site_pattern",
      "participant_day",
      run_id
    )
  )
  site_test <- h02_joint_smooth_wald(
    candidates$site_pattern,
    "s(time_hour,site)"
  )
  site_row <- comparisons$comparison_id == "site_pattern_vs_no_site"
  comparisons$test_method[site_row] <- paste(
    "joint approximate Bayesian Wald chi-square test of all sum-to-zero",
    "site-smooth coefficients"
  )
  comparisons$test_statistic[site_row] <- site_test["statistic"]
  comparisons$test_df[site_row] <- site_test["df"]
  comparisons$p_raw[site_row] <- site_test["p"]
  comparisons$log_p_raw[site_row] <- site_test["log_p"]
  selected_id <- "site_pattern"
  final <- candidates[[selected_id]]
  list(
    data = data,
    preliminary = preliminary,
    candidates = candidates,
    final = final,
    rho = rho,
    selected_model_id = selected_id,
    model_table = model_table,
    comparisons = comparisons
  )
}

h02_fit_selected_run <- function(
  frame,
  run_id,
  selected_model_id
) {
  data <- h02_prepare_fit_data(frame)
  formula <- h02_formula_set()[[selected_model_id]]
  preliminary <- h02_fit_bam(
    formula,
    data,
    method = "fREML",
    rho = 0
  )
  rho <- h02_estimate_rho(preliminary, data)
  final <- h02_fit_bam(
    formula,
    data,
    method = "fREML",
    rho = rho
  )
  list(
    data = data,
    preliminary = preliminary,
    final = final,
    rho = rho,
    selected_model_id = selected_model_id,
    model_table = h02_model_row(
      final,
      selected_model_id,
      run_id,
      rho
    )
  )
}

h02_residual_acf <- function(fit, data, stage, run_id, max_lag = 6L) {
  residual <- if (
    startsWith(stage, "final") &&
      !is.null(fit$std.rsd) &&
      length(fit$std.rsd) == nrow(data)
  ) {
    fit$std.rsd
  } else {
    stats::residuals(fit, type = "response")
  }
  dplyr::bind_rows(lapply(seq_len(max_lag), function(lag) {
    estimate <- h02_boundary_lag_correlation(
      residual,
      data$AR_start,
      lag
    )
    tibble::tibble(
      run_id = run_id,
      stage = stage,
      lag_30_minute_bins = lag,
      correlation = unname(estimate["correlation"]),
      pairs = as.integer(estimate["pairs"])
    )
  }))
}

h02_residual_summary <- function(fit, data, run_id) {
  residual <- if (
    !is.null(fit$std.rsd) &&
      length(fit$std.rsd) == nrow(data)
  ) {
    fit$std.rsd
  } else {
    stats::residuals(fit, type = "response")
  }
  fitted <- stats::fitted(fit)
  q <- stats::quantile(
    residual,
    probs = c(0.01, 0.05, 0.5, 0.95, 0.99),
    na.rm = TRUE,
    names = FALSE
  )
  tibble::tibble(
    run_id = run_id,
    residual_type = if (!is.null(fit$std.rsd)) {
      "AR-standardized"
    } else {
      "response"
    },
    n = length(residual),
    mean = mean(residual),
    sd = stats::sd(residual),
    rmse = sqrt(mean(residual^2)),
    q01 = q[1],
    q05 = q[2],
    median = q[3],
    q95 = q[4],
    q99 = q[5],
    correlation_absolute_residual_fitted = stats::cor(
      abs(residual),
      fitted,
      use = "complete.obs"
    ),
    maximum_absolute_residual = max(abs(residual), na.rm = TRUE)
  )
}

h02_boundary_audit <- function(data, run_id) {
  data |>
    dplyr::count(.data$ar_start_reason, name = "observations") |>
    dplyr::mutate(
      run_id = run_id,
      AR_start = .data$ar_start_reason != "continuous",
      .before = 1L
    )
}

h02_smooth_labels <- function(fit) {
  vapply(fit$smooth, function(smooth) smooth$label, character(1))
}

h02_random_smooth_indices <- function(fit) {
  labels <- h02_smooth_labels(fit)
  random <- grepl("participant", labels, fixed = TRUE)
  unique(unlist(lapply(fit$smooth[random], function(smooth) {
    seq.int(smooth$first.para, smooth$last.para)
  })))
}

h02_find_smooth <- function(fit, pattern) {
  labels <- h02_smooth_labels(fit)
  which(grepl(pattern, labels, fixed = TRUE))
}

h02_site_predictions <- function(
  fit,
  data,
  run_id,
  simulations = 2000L,
  seed = h02_seed(run_id, 10L)
) {
  sites <- levels(data$site)
  time_hour <- (seq.int(0L, 1410L, by = 30L) + 15) / 60
  grid <- tidyr::crossing(
    site = factor(sites, levels = sites),
    time_hour = time_hour
  ) |>
    dplyr::mutate(
      site_smooth = ordered(as.character(.data$site), levels = sites),
      participant = factor(
        levels(data$participant)[1],
        levels = levels(data$participant)
      ),
      participant_day = factor(
        levels(data$participant_day)[1],
        levels = levels(data$participant_day)
      )
    )

  # predict.bam()'s discrete new-data path fails for `sz` factor smooths
  # (mgcv 1.9-4) when constructing an lpmatrix. The parent GAM method uses
  # the same fitted coefficients and basis definitions without that failure.
  L <- mgcv::predict.gam(
    fit,
    newdata = as.data.frame(grid),
    type = "lpmatrix",
    newdata.guaranteed = TRUE
  )
  random_columns <- h02_random_smooth_indices(fit)
  if (length(random_columns) > 0L) {
    L[, random_columns] <- 0
  }
  keep <- colSums(abs(L)) > 0
  L_keep <- L[, keep, drop = FALSE]
  coefficients <- stats::coef(fit)[keep]
  covariance <- fit$Vp[keep, keep, drop = FALSE]
  eta <- drop(L_keep %*% coefficients)
  se <- sqrt(pmax(0, rowSums((L_keep %*% covariance) * L_keep)))

  set.seed(seed)
  draws <- mgcv::rmvn(
    simulations,
    mu = rep(0, length(coefficients)),
    V = covariance
  )
  prediction_error <- L_keep %*% t(draws)
  standardized <- sweep(
    prediction_error,
    1L,
    pmax(se, sqrt(.Machine$double.eps)),
    "/"
  )
  critical <- unname(stats::quantile(
    apply(abs(standardized), 2L, max),
    probs = 0.95,
    names = FALSE
  ))
  grid <- grid |>
    dplyr::mutate(
      run_id = run_id,
      eta = eta,
      standard_error = se,
      simultaneous_critical_value = critical,
      lower_eta = .data$eta - critical * .data$standard_error,
      upper_eta = .data$eta + critical * .data$standard_error,
      estimate_melEDI_lx = h02_inverse_transform(.data$eta),
      lower_melEDI_lx = h02_inverse_transform(.data$lower_eta),
      upper_melEDI_lx = h02_inverse_transform(.data$upper_eta),
      clock_bin = as.integer(round(.data$time_hour * 60 - 15))
    )

  mean_rows <- split(seq_len(nrow(grid)), grid$time_hour)
  L_mean <- do.call(
    rbind,
    lapply(mean_rows, function(rows) {
      colMeans(L_keep[rows, , drop = FALSE])
    })
  )
  mean_lookup <- match(grid$time_hour, sort(unique(grid$time_hour)))
  L_deviation <- L_keep - L_mean[mean_lookup, , drop = FALSE]
  deviation_eta <- drop(L_deviation %*% coefficients)
  deviation_se <- sqrt(pmax(
    0,
    rowSums((L_deviation %*% covariance) * L_deviation)
  ))
  deviation_error <- L_deviation %*% t(draws)
  deviation_standardized <- sweep(
    deviation_error,
    1L,
    pmax(deviation_se, sqrt(.Machine$double.eps)),
    "/"
  )
  deviation_critical <- unname(stats::quantile(
    apply(abs(deviation_standardized), 2L, max),
    probs = 0.95,
    names = FALSE
  ))
  grid |>
    dplyr::mutate(
      equal_site_deviation_eta = deviation_eta,
      deviation_standard_error = deviation_se,
      deviation_simultaneous_critical_value = deviation_critical,
      site_to_equal_site_ratio = 10^.data$equal_site_deviation_eta,
      ratio_lower = 10^(.data$equal_site_deviation_eta -
        deviation_critical * .data$deviation_standard_error),
      ratio_upper = 10^(.data$equal_site_deviation_eta +
        deviation_critical * .data$deviation_standard_error)
    ) |>
    dplyr::select(
      "run_id",
      "site",
      "clock_bin",
      "time_hour",
      "eta",
      "standard_error",
      "lower_eta",
      "upper_eta",
      "estimate_melEDI_lx",
      "lower_melEDI_lx",
      "upper_melEDI_lx",
      "simultaneous_critical_value",
      "equal_site_deviation_eta",
      "deviation_standard_error",
      "site_to_equal_site_ratio",
      "ratio_lower",
      "ratio_upper",
      "deviation_simultaneous_critical_value"
    )
}

h02_site_windows <- function(site_predictions, run_id) {
  selected <- site_predictions |>
    dplyr::filter(.data$run_id == .env$run_id) |>
    dplyr::arrange(.data$site, .data$clock_bin) |>
    dplyr::mutate(
      direction = dplyr::case_when(
        .data$ratio_lower > 1 ~ "higher",
        .data$ratio_upper < 1 ~ "lower",
        TRUE ~ "no_simultaneous_difference"
      )
    ) |>
    dplyr::group_by(.data$site) |>
    dplyr::mutate(
      window_start = dplyr::row_number() == 1L |
        .data$direction !=
          dplyr::lag(
            .data$direction,
            default = dplyr::first(.data$direction)
          ),
      window_id = cumsum(.data$window_start)
    ) |>
    dplyr::ungroup()

  windows <- selected |>
    dplyr::filter(.data$direction != "no_simultaneous_difference") |>
    dplyr::group_by(.data$site, .data$window_id, .data$direction) |>
    dplyr::summarise(
      start_clock_bin = min(.data$clock_bin),
      end_clock_bin_exclusive = max(.data$clock_bin) + 30L,
      bins_30_minute = dplyr::n(),
      minimum_point_ratio = min(.data$site_to_equal_site_ratio),
      maximum_point_ratio = max(.data$site_to_equal_site_ratio),
      minimum_simultaneous_lower = min(.data$ratio_lower),
      maximum_simultaneous_upper = max(.data$ratio_upper),
      .groups = "drop"
    )
  no_windows <- setdiff(
    unique(as.character(selected$site)),
    unique(as.character(windows$site))
  )
  if (length(no_windows) > 0L) {
    windows <- dplyr::bind_rows(
      windows,
      tibble::tibble(
        site = no_windows,
        window_id = NA_integer_,
        direction = "no_simultaneous_difference",
        start_clock_bin = NA_integer_,
        end_clock_bin_exclusive = NA_integer_,
        bins_30_minute = 0L,
        minimum_point_ratio = NA_real_,
        maximum_point_ratio = NA_real_,
        minimum_simultaneous_lower = NA_real_,
        maximum_simultaneous_upper = NA_real_
      )
    )
  }
  format_clock <- function(minutes) {
    ifelse(
      is.na(minutes),
      NA_character_,
      sprintf("%02d:%02d", (minutes %/% 60L) %% 24L, minutes %% 60L)
    )
  }
  windows |>
    dplyr::mutate(
      run_id = run_id,
      start_local_clock = format_clock(.data$start_clock_bin),
      end_local_clock = dplyr::if_else(
        .data$end_clock_bin_exclusive == 1440L,
        "24:00",
        format_clock(.data$end_clock_bin_exclusive)
      ),
      interval_definition = paste(
        "contiguous 30-minute bins whose simultaneous 95% interval for",
        "the site-to-equal-site ratio excludes 1"
      ),
      .before = 1L
    ) |>
    dplyr::arrange(.data$site, .data$start_clock_bin)
}

h02_term_contribution <- function(fit, newdata, smooth_index) {
  smooth <- fit$smooth[[smooth_index]]
  columns <- seq.int(smooth$first.para, smooth$last.para)
  L <- mgcv::predict.gam(
    fit,
    newdata = as.data.frame(newdata),
    type = "lpmatrix",
    newdata.guaranteed = TRUE
  )
  drop(L[, columns, drop = FALSE] %*% stats::coef(fit)[columns])
}

h02_fitted_contributions <- function(fit, data, site_predictions) {
  participant_index <- h02_find_smooth(
    fit,
    "s(time_hour,participant)"
  )
  day_index <- h02_find_smooth(fit, "s(participant_day)")
  if (length(participant_index) != 1L || length(day_index) != 1L) {
    h02_abort(
      "Could not uniquely identify participant and participant-day smooths"
    )
  }
  participant_info <- data |>
    dplyr::distinct(.data$site, .data$participant) |>
    dplyr::arrange(.data$site, .data$participant)
  participant_grid <- tidyr::crossing(
    participant_info,
    time_hour = (seq.int(0L, 1410L, by = 30L) + 15) / 60
  ) |>
    dplyr::mutate(
      site_smooth = ordered(
        as.character(.data$site),
        levels = levels(data$site)
      ),
      participant_day = factor(
        levels(data$participant_day)[1],
        levels = levels(data$participant_day)
      )
    )
  participant_grid$participant_eta <- h02_term_contribution(
    fit,
    participant_grid,
    participant_index
  )

  day_grid <- data |>
    dplyr::distinct(
      .data$site,
      .data$participant,
      .data$participant_day
    ) |>
    dplyr::mutate(
      site_smooth = ordered(
        as.character(.data$site),
        levels = levels(data$site)
      ),
      time_hour = 0.25
    )
  day_grid$participant_day_eta <- h02_term_contribution(
    fit,
    day_grid,
    day_index
  )
  list(
    site = site_predictions |>
      dplyr::select("site", "time_hour", "eta"),
    participant = participant_grid |>
      dplyr::transmute(
        site = as.character(.data$site),
        participant = as.character(.data$participant),
        time_hour = .data$time_hour,
        participant_eta = .data$participant_eta
      ),
    participant_day = day_grid |>
      dplyr::transmute(
        site = as.character(.data$site),
        participant = as.character(.data$participant),
        participant_day = as.character(.data$participant_day),
        participant_day_eta = .data$participant_day_eta
      )
  )
}

h02_variation_point <- function(contributions) {
  site_value <- contributions$site |>
    dplyr::group_by(.data$time_hour) |>
    dplyr::summarise(value = stats::var(.data$eta), .groups = "drop") |>
    dplyr::summarise(value = mean(.data$value)) |>
    dplyr::pull(.data$value)

  participant_value <- contributions$participant |>
    dplyr::group_by(.data$site, .data$time_hour) |>
    dplyr::summarise(
      value = stats::var(.data$participant_eta),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$site) |>
    dplyr::summarise(value = mean(.data$value), .groups = "drop") |>
    dplyr::summarise(value = mean(.data$value)) |>
    dplyr::pull(.data$value)

  day_value <- contributions$participant_day |>
    dplyr::group_by(.data$site, .data$participant) |>
    dplyr::summarise(
      value = if (dplyr::n() > 1L) {
        stats::var(.data$participant_day_eta)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$site) |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::summarise(value = mean(.data$value)) |>
    dplyr::pull(.data$value)

  c(
    site_curve_variation = site_value,
    participant_curve_variation = participant_value,
    participant_day_intercept_variation = day_value,
    participant_plus_day_variation = participant_value + day_value,
    participant_to_site_ratio = participant_value / site_value,
    participant_plus_day_to_site_ratio = (participant_value + day_value) /
      site_value
  )
}

h02_bootstrap_variation <- function(
  contributions,
  run_id,
  replicates = 2000L,
  seed = h02_seed(run_id, 20L)
) {
  site_names <- sort(unique(as.character(contributions$site$site)))
  time_values <- sort(unique(contributions$site$time_hour))
  site_matrix <- tidyr::pivot_wider(
    contributions$site,
    names_from = "time_hour",
    values_from = "eta"
  )
  site_matrix <- as.matrix(
    site_matrix[
      match(site_names, as.character(site_matrix$site)),
      -1,
      drop = FALSE
    ]
  )

  participant_by_site <- lapply(site_names, function(site_name) {
    x <- dplyr::filter(
      contributions$participant,
      .data$site == site_name
    ) |>
      tidyr::pivot_wider(
        names_from = "time_hour",
        values_from = "participant_eta"
      ) |>
      dplyr::arrange(.data$participant)
    matrix <- as.matrix(x[, as.character(time_values), drop = FALSE])
    rownames(matrix) <- x$participant
    matrix
  })
  names(participant_by_site) <- site_names

  day_by_site <- lapply(site_names, function(site_name) {
    x <- dplyr::filter(
      contributions$participant_day,
      .data$site == site_name
    )
    split(x$participant_day_eta, x$participant)
  })
  names(day_by_site) <- site_names

  set.seed(seed)
  output <- matrix(
    NA_real_,
    nrow = replicates,
    ncol = 6L,
    dimnames = list(
      NULL,
      c(
        "site_curve_variation",
        "participant_curve_variation",
        "participant_day_intercept_variation",
        "participant_plus_day_variation",
        "participant_to_site_ratio",
        "participant_plus_day_to_site_ratio"
      )
    )
  )
  for (b in seq_len(replicates)) {
    sampled_sites <- sample(
      seq_along(site_names),
      length(site_names),
      replace = TRUE
    )
    site_value <- mean(apply(
      site_matrix[sampled_sites, , drop = FALSE],
      2L,
      stats::var
    ))
    participant_site_values <- numeric(length(sampled_sites))
    day_site_values <- numeric(length(sampled_sites))
    for (j in seq_along(sampled_sites)) {
      site_index <- sampled_sites[j]
      participant_matrix <- participant_by_site[[site_index]]
      participant_ids <- rownames(participant_matrix)
      sampled_participant_rows <- sample(
        seq_along(participant_ids),
        length(participant_ids),
        replace = TRUE
      )
      participant_site_values[j] <- mean(apply(
        participant_matrix[
          sampled_participant_rows,
          ,
          drop = FALSE
        ],
        2L,
        stats::var
      ))
      sampled_participant_ids <- participant_ids[
        sampled_participant_rows
      ]
      day_variances <- vapply(
        sampled_participant_ids,
        function(participant_id) {
          values <- day_by_site[[site_index]][[participant_id]]
          if (length(values) < 2L) {
            return(NA_real_)
          }
          stats::var(sample(values, length(values), replace = TRUE))
        },
        numeric(1)
      )
      day_site_values[j] <- mean(day_variances, na.rm = TRUE)
    }
    participant_value <- mean(participant_site_values)
    day_value <- mean(day_site_values)
    output[b, ] <- c(
      site_value,
      participant_value,
      day_value,
      participant_value + day_value,
      participant_value / site_value,
      (participant_value + day_value) / site_value
    )
  }
  output
}

h02_variation_summary <- function(
  fit,
  data,
  site_predictions,
  run_id,
  replicates = bootstrap_count(2000L)
) {
  contributions <- h02_fitted_contributions(
    fit,
    data,
    site_predictions
  )
  point <- h02_variation_point(contributions)
  bootstrap <- h02_bootstrap_variation(
    contributions,
    run_id,
    replicates = replicates
  )
  intervals <- apply(
    bootstrap,
    2L,
    stats::quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  definitions <- c(
    site_curve_variation = paste(
      "Mean across 48 clock bins of the sample variance across equal-weight",
      "site fitted curves"
    ),
    participant_curve_variation = paste(
      "Equal-site mean across 48 clock bins of within-site sample variance",
      "among fitted participant deviation curves"
    ),
    participant_day_intercept_variation = paste(
      "Equal-site mean of participant-specific sample variance among fitted",
      "participant-day intercept contributions"
    ),
    participant_plus_day_variation = paste(
      "Sum of fitted participant-curve and participant-day-intercept",
      "variation summaries"
    ),
    participant_to_site_ratio = "Participant-curve variation divided by site-curve variation",
    participant_plus_day_to_site_ratio = paste(
      "Participant-curve plus participant-day-intercept variation divided",
      "by site-curve variation"
    )
  )
  units <- c(
    rep("squared log10(melEDI + 0.1 lx) prediction units", 4L),
    "ratio",
    "ratio"
  )
  summary <- tibble::tibble(
    run_id = .env$run_id,
    summary_id = names(point),
    definition = unname(definitions[names(point)]),
    estimate = unname(point),
    lower_95 = intervals[1L, names(point)],
    upper_95 = intervals[2L, names(point)],
    unit = unname(units),
    confidence_interval_method = paste(
      "percentile hierarchical cluster bootstrap of fitted contributions;",
      "sites, participants within sites, and days within participants;",
      replicates,
      "replicates; conditional on fitted smoothing structure"
    ),
    bootstrap_seed = h02_seed(.env$run_id, 20L),
    bootstrap_replicates = replicates
  )
  list(
    summary = summary,
    contributions = contributions,
    bootstrap = bootstrap
  )
}

h02_contribution_influence <- function(
  contributions,
  reference,
  run_id
) {
  reference_named <- stats::setNames(
    reference$estimate,
    reference$summary_id
  )
  sites <- sort(unique(as.character(contributions$site$site)))
  site_results <- dplyr::bind_rows(lapply(sites, function(site_name) {
    reduced <- list(
      site = dplyr::filter(
        contributions$site,
        as.character(.data$site) != site_name
      ),
      participant = dplyr::filter(
        contributions$participant,
        .data$site != site_name
      ),
      participant_day = dplyr::filter(
        contributions$participant_day,
        .data$site != site_name
      )
    )
    value <- h02_variation_point(reduced)
    tibble::tibble(
      deleted_level = "site",
      deleted_site = site_name,
      deleted_participant = NA_character_,
      summary_id = names(value),
      deletion_estimate = unname(value)
    )
  }))

  participants <- contributions$participant |>
    dplyr::distinct(.data$site, .data$participant)
  participant_results <- dplyr::bind_rows(lapply(
    seq_len(nrow(participants)),
    function(i) {
      participant_id <- participants$participant[i]
      site_name <- participants$site[i]
      reduced <- list(
        site = contributions$site,
        participant = dplyr::filter(
          contributions$participant,
          .data$participant != participant_id
        ),
        participant_day = dplyr::filter(
          contributions$participant_day,
          .data$participant != participant_id
        )
      )
      value <- h02_variation_point(reduced)
      tibble::tibble(
        deleted_level = "participant",
        deleted_site = site_name,
        deleted_participant = participant_id,
        summary_id = names(value),
        deletion_estimate = unname(value)
      )
    }
  ))
  dplyr::bind_rows(site_results, participant_results) |>
    dplyr::mutate(
      run_id = run_id,
      reference_estimate = unname(reference_named[.data$summary_id]),
      relative_change = dplyr::if_else(
        .data$reference_estimate == 0,
        NA_real_,
        (.data$deletion_estimate - .data$reference_estimate) /
          .data$reference_estimate
      ),
      diagnostic_scope = paste(
        "conditional deletion from fitted contributions;",
        "the GAMM was not refitted"
      ),
      .before = 1L
    )
}

h02_k_check <- function(fit, run_id) {
  diagnostic_seed <- h02_seed(run_id, 30L)
  set.seed(diagnostic_seed)
  result <- mgcv::k.check(fit, subsample = 5000L, n.rep = 200L)
  tibble::as_tibble(result, rownames = "smooth") |>
    dplyr::rename(
      k_prime = "k'",
      edf = "edf",
      k_index = "k-index",
      p_value = "p-value"
    ) |>
    dplyr::mutate(
      run_id = run_id,
      diagnostic_seed = diagnostic_seed,
      simulation_replicates = 200L,
      .before = 1L
    )
}

h02_influence_scores <- function(fit, data, run_id) {
  residual <- if (!is.null(fit$std.rsd)) {
    fit$std.rsd
  } else {
    stats::residuals(fit, type = "response")
  }
  leverage <- if (!is.null(fit$hat)) {
    fit$hat
  } else {
    rep(NA_real_, nrow(data))
  }
  scored <- data |>
    dplyr::transmute(
      site = as.character(.data$site),
      participant = as.character(.data$participant),
      residual = residual,
      leverage = leverage
    )
  participant <- scored |>
    dplyr::group_by(.data$site, .data$participant) |>
    dplyr::summarise(
      observations = dplyr::n(),
      mean_absolute_standardized_residual = mean(abs(.data$residual)),
      p99_absolute_standardized_residual = unname(stats::quantile(
        abs(.data$residual),
        0.99
      )),
      maximum_absolute_standardized_residual = max(abs(.data$residual)),
      total_leverage = if (any(is.finite(.data$leverage))) {
        sum(.data$leverage, na.rm = TRUE)
      } else {
        NA_real_
      },
      maximum_leverage = if (any(is.finite(.data$leverage))) {
        max(.data$leverage, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      run_id = run_id,
      level = "participant",
      .before = 1L
    )
  site <- scored |>
    dplyr::group_by(.data$site) |>
    dplyr::summarise(
      observations = dplyr::n(),
      mean_absolute_standardized_residual = mean(abs(.data$residual)),
      p99_absolute_standardized_residual = unname(stats::quantile(
        abs(.data$residual),
        0.99
      )),
      maximum_absolute_standardized_residual = max(abs(.data$residual)),
      total_leverage = if (any(is.finite(.data$leverage))) {
        sum(.data$leverage, na.rm = TRUE)
      } else {
        NA_real_
      },
      maximum_leverage = if (any(is.finite(.data$leverage))) {
        max(.data$leverage, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      run_id = run_id,
      level = "site",
      participant = NA_character_,
      .before = 1L
    )
  dplyr::bind_rows(participant, site)
}

h02_selected_model_specification <- function(
  fit,
  selected_model_id,
  rho,
  run_id
) {
  spec <- h02_specification_table()
  extra <- tibble::tribble(
    ~field,
    ~value,
    "selected_model_id",
    selected_model_id,
    "selected_formula",
    paste(deparse(stats::formula(fit)), collapse = " "),
    "selected_primary_rho",
    format(rho, digits = 17),
    "selected_from_run_id",
    run_id,
    "mgcv_version",
    as.character(utils::packageVersion("mgcv")),
    "R_version",
    as.character(getRversion()),
    "inheritance_rule_for_H11",
    paste(
      "Use this formula, transform, basis dimensions, knots, estimation",
      "method, and AR-boundary algorithm verbatim; re-estimate rho from",
      "boundary-aware preliminary residuals in the inherited model frame."
    )
  )
  dplyr::bind_rows(spec, extra)
}
