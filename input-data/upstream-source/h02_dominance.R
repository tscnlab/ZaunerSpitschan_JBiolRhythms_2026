# H02 exact conditional Shapley/general-dominance allocation.

h02_dominance_registry <- function() {
  tibble::tribble(
    ~run_id,
    ~placement,
    ~analytical_role,
    "main__glasses__all_available",
    "glasses",
    "primary_near_eye",
    "main__chest__all_available",
    "chest",
    "complementary_chest"
  )
}

h02_dominance_base_term <- function(spec = h02_model_specification()) {
  sprintf("s(time_hour, bs = 'cc', k = %d)", spec$overall_k)
}

h02_dominance_terms <- function(spec = h02_model_specification()) {
  c(
    site_pattern = sprintf(
      "s(time_hour, site, bs = 'sz', k = %d)",
      spec$site_pattern_k
    ),
    participant_pattern = sprintf(
      "s(time_hour, participant, bs = 'fs', k = %d)",
      spec$participant_k
    ),
    participant_day = "s(participant_day, bs = 're')"
  )
}

h02_dominance_design <- function(
  base_term = h02_dominance_base_term(),
  terms = h02_dominance_terms()
) {
  component_names <- names(terms)
  if (
    length(component_names) == 0L ||
      any(!nzchar(component_names)) ||
      anyDuplicated(component_names)
  ) {
    h02_abort("Dominance components must have unique non-empty names")
  }
  component_count <- length(component_names)
  masks <- seq.int(0L, bitwShiftL(1L, component_count) - 1L)
  included <- lapply(masks, function(mask) {
    component_names[vapply(
      seq_along(component_names),
      function(index) {
        bitwAnd(mask, bitwShiftL(1L, index - 1L)) != 0L
      },
      logical(1)
    )]
  })
  formulas <- lapply(included, function(components) {
    rhs <- c(base_term, unname(terms[components]))
    stats::as.formula(
      paste("response ~", paste(rhs, collapse = " + "))
    )
  })
  tibble::tibble(
    mask = masks,
    subset_size = lengths(included),
    subset_id = vapply(
      included,
      function(x) {
        if (length(x) == 0L) {
          "common_time_only"
        } else {
          paste(c("common_time", x), collapse = "+")
        }
      },
      character(1)
    ),
    included_components = vapply(
      included,
      paste,
      collapse = ";",
      character(1)
    ),
    formula = formulas,
    formula_text = vapply(
      formulas,
      function(x) paste(deparse(x), collapse = " "),
      character(1)
    )
  )
}

h02_dominance_map <- function(design, component_names) {
  component_count <- length(component_names)
  expected_masks <- seq.int(
    0L,
    bitwShiftL(1L, component_count) - 1L
  )
  if (!identical(sort(design$mask), expected_masks)) {
    h02_abort("Dominance design does not contain every required subset")
  }
  dplyr::bind_rows(lapply(seq_along(component_names), function(index) {
    bit <- bitwShiftL(1L, index - 1L)
    without <- design |>
      dplyr::filter(bitwAnd(.data$mask, bit) == 0L)
    with_mask <- bitwOr(without$mask, bit)
    with_index <- match(with_mask, design$mask)
    if (anyNA(with_index)) {
      h02_abort("Dominance design is missing required component subsets")
    }
    subset_size <- without$subset_size
    tibble::tibble(
      component = component_names[index],
      without_mask = without$mask,
      with_mask = with_mask,
      without_subset_id = without$subset_id,
      with_subset_id = design$subset_id[with_index],
      subset_size = subset_size,
      shapley_weight = factorial(subset_size) *
        factorial(component_count - subset_size - 1L) /
        factorial(component_count)
    )
  }))
}

h02_dominance_values <- function(response, predictions, design) {
  if (nrow(predictions) != length(response)) {
    h02_abort("Dominance predictions and response have incompatible rows")
  }
  if (ncol(predictions) != nrow(design)) {
    h02_abort(
      "Dominance predictions and subset design have incompatible columns"
    )
  }
  if (any(!is.finite(response)) || any(!is.finite(predictions))) {
    h02_abort("Dominance response and predictions must be finite")
  }
  centered_response <- response - mean(response)
  total_sum_squares <- sum(centered_response^2)
  if (!is.finite(total_sum_squares) || total_sum_squares <= 0) {
    h02_abort("Dominance response total sum of squares is not positive")
  }
  residual <- sweep(predictions, 1L, response, "-")
  squared_error <- colSums(residual^2)
  names(squared_error) <- design$mask
  r_squared <- 1 - squared_error / total_sum_squares
  names(r_squared) <- design$mask
  list(
    r_squared = r_squared,
    squared_error = squared_error,
    total_sum_squares = total_sum_squares
  )
}

h02_dominance_comparisons <- function(
  heterogeneity_summary,
  heterogeneity_increment_R2
) {
  allocation <- stats::setNames(
    heterogeneity_summary$allocated_R2,
    heterogeneity_summary$component
  )
  site <- unname(allocation["site_pattern"])
  participant <- unname(allocation["participant_pattern"])
  day <- unname(allocation["participant_day"])
  estimates <- c(
    participant_to_site_shapley_ratio = participant / site,
    participant_plus_day_to_site_shapley_ratio = (participant + day) / site,
    participant_share_of_site_plus_participant = participant /
      (site + participant),
    participant_plus_day_share_of_heterogeneity = (participant + day) /
      heterogeneity_increment_R2
  )
  definitions <- c(
    participant_to_site_shapley_ratio = paste(
      "Conditional Shapley R2 allocated to the participant-pattern block",
      "divided by that allocated to the site-pattern block"
    ),
    participant_plus_day_to_site_shapley_ratio = paste(
      "Conditional Shapley R2 allocated to participant-pattern plus",
      "participant-day blocks divided by that allocated to the",
      "site-pattern block"
    ),
    participant_share_of_site_plus_participant = paste(
      "Participant-pattern allocation divided by the sum of site-pattern",
      "and participant-pattern allocations"
    ),
    participant_plus_day_share_of_heterogeneity = paste(
      "Participant-pattern plus participant-day allocation divided by",
      "the full-model R2 increment beyond the common time curve"
    )
  )
  tibble::tibble(
    comparison_id = names(estimates),
    definition = unname(definitions[names(estimates)]),
    estimate = unname(estimates),
    unit = c("ratio", "ratio", "proportion", "proportion")
  )
}

h02_dominance_decomposition <- function(values, design, dominance_map) {
  r_squared <- values$r_squared
  without_index <- match(dominance_map$without_mask, design$mask)
  with_index <- match(dominance_map$with_mask, design$mask)
  marginal <- dominance_map |>
    dplyr::mutate(
      R2_without = unname(r_squared[without_index]),
      R2_with = unname(r_squared[with_index]),
      incremental_R2 = .data$R2_with - .data$R2_without,
      weighted_increment = .data$shapley_weight * .data$incremental_R2
    )
  full_mask <- max(design$mask)
  common_time_R2 <- unname(r_squared[match(0L, design$mask)])
  full_model_R2 <- unname(r_squared[match(full_mask, design$mask)])
  heterogeneity_increment_R2 <- full_model_R2 - common_time_R2
  heterogeneity_summary <- marginal |>
    dplyr::group_by(.data$component) |>
    dplyr::summarise(
      allocated_R2 = sum(.data$weighted_increment),
      minimum_incremental_R2 = min(.data$incremental_R2),
      maximum_incremental_R2 = max(.data$incremental_R2),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      allocation_method = "exact conditional Shapley/general dominance beyond common time",
      general_dominance_rank = rank(
        -.data$allocated_R2,
        ties.method = "min"
      )
    )
  tolerance <- 100 *
    .Machine$double.eps *
    max(
      1,
      abs(common_time_R2),
      abs(full_model_R2)
    )
  if (
    !isTRUE(all.equal(
      sum(heterogeneity_summary$allocated_R2),
      heterogeneity_increment_R2,
      tolerance = tolerance
    ))
  ) {
    h02_abort(
      "Conditional Shapley allocations do not sum to the full-model R2 increment"
    )
  }
  allocation <- dplyr::bind_rows(
    tibble::tibble(
      component = "common_time",
      allocated_R2 = common_time_R2,
      minimum_incremental_R2 = NA_real_,
      maximum_incremental_R2 = NA_real_,
      allocation_method = "mandatory common-time baseline fitted before heterogeneity blocks",
      general_dominance_rank = NA_real_
    ),
    heterogeneity_summary
  ) |>
    dplyr::mutate(
      share_of_full_model_R2 = .data$allocated_R2 / full_model_R2,
      share_of_increment_beyond_common = dplyr::if_else(
        .data$component == "common_time",
        NA_real_,
        .data$allocated_R2 / heterogeneity_increment_R2
      ),
      common_time_R2 = common_time_R2,
      heterogeneity_increment_R2 = heterogeneity_increment_R2,
      full_model_R2 = full_model_R2
    )
  conditional <- marginal |>
    dplyr::group_by(.data$component, .data$subset_size) |>
    dplyr::summarise(
      conditional_dominance_R2 = mean(.data$incremental_R2),
      subsets = dplyr::n(),
      .groups = "drop"
    )
  comparisons <- h02_dominance_comparisons(
    heterogeneity_summary,
    heterogeneity_increment_R2
  )
  list(
    marginal = marginal,
    conditional = conditional,
    allocation = allocation,
    comparisons = comparisons
  )
}

h02_dominance_hierarchy <- function(data) {
  required <- c("site", "participant", "participant_day")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    h02_abort(
      "Dominance hierarchy is missing required columns: %s",
      paste(missing, collapse = ", ")
    )
  }
  site_names <- sort(unique(as.character(data$site)))
  hierarchy <- lapply(site_names, function(site_name) {
    site_rows <- which(as.character(data$site) == site_name)
    participant_names <- sort(unique(as.character(
      data$participant[site_rows]
    )))
    participants <- lapply(participant_names, function(participant_name) {
      participant_rows <- site_rows[
        as.character(data$participant[site_rows]) == participant_name
      ]
      day_names <- sort(unique(as.character(
        data$participant_day[participant_rows]
      )))
      days <- lapply(day_names, function(day_name) {
        participant_rows[
          as.character(data$participant_day[participant_rows]) == day_name
        ]
      })
      names(days) <- day_names
      days
    })
    names(participants) <- participant_names
    participants
  })
  names(hierarchy) <- site_names
  hierarchy
}

h02_sample_dominance_rows <- function(hierarchy) {
  if (length(hierarchy) == 0L) {
    h02_abort("Dominance hierarchy contains no sites")
  }
  sampled_site_positions <- sample.int(
    length(hierarchy),
    length(hierarchy),
    replace = TRUE
  )
  unlist(
    lapply(sampled_site_positions, function(site_position) {
      participants <- hierarchy[[site_position]]
      sampled_participant_positions <- sample.int(
        length(participants),
        length(participants),
        replace = TRUE
      )
      unlist(
        lapply(
          sampled_participant_positions,
          function(participant_position) {
            days <- participants[[participant_position]]
            sampled_day_positions <- sample.int(
              length(days),
              length(days),
              replace = TRUE
            )
            unlist(days[sampled_day_positions], use.names = FALSE)
          }
        ),
        use.names = FALSE
      )
    }),
    use.names = FALSE
  )
}

h02_bootstrap_dominance <- function(
  response,
  predictions,
  data,
  design,
  dominance_map,
  replicates = 2000L,
  seed
) {
  if (
    length(replicates) != 1L ||
      is.na(replicates) ||
      replicates < 1L ||
      replicates != as.integer(replicates)
  ) {
    h02_abort("Dominance bootstrap replicates must be a positive integer")
  }
  component_names <- c(
    "common_time",
    unique(dominance_map$component)
  )
  comparison_names <- c(
    "participant_to_site_shapley_ratio",
    "participant_plus_day_to_site_shapley_ratio",
    "participant_share_of_site_plus_participant",
    "participant_plus_day_share_of_heterogeneity"
  )
  hierarchy <- h02_dominance_hierarchy(data)
  allocated_R2 <- matrix(
    NA_real_,
    nrow = replicates,
    ncol = length(component_names),
    dimnames = list(NULL, component_names)
  )
  share_full <- allocated_R2
  share_increment <- matrix(
    NA_real_,
    nrow = replicates,
    ncol = length(component_names) - 1L,
    dimnames = list(NULL, component_names[-1L])
  )
  comparison <- matrix(
    NA_real_,
    nrow = replicates,
    ncol = length(comparison_names),
    dimnames = list(NULL, comparison_names)
  )
  common_time_R2 <- rep(NA_real_, replicates)
  heterogeneity_increment_R2 <- rep(NA_real_, replicates)
  full_model_R2 <- rep(NA_real_, replicates)
  set.seed(seed)
  for (replicate in seq_len(replicates)) {
    rows <- h02_sample_dominance_rows(hierarchy)
    values <- h02_dominance_values(
      response[rows],
      predictions[rows, , drop = FALSE],
      design
    )
    decomposition <- h02_dominance_decomposition(
      values,
      design,
      dominance_map
    )
    allocation_index <- match(
      component_names,
      decomposition$allocation$component
    )
    comparison_index <- match(
      comparison_names,
      decomposition$comparisons$comparison_id
    )
    allocated_R2[replicate, ] <-
      decomposition$allocation$allocated_R2[allocation_index]
    share_full[replicate, ] <-
      decomposition$allocation$share_of_full_model_R2[allocation_index]
    share_increment[replicate, ] <-
      decomposition$allocation$share_of_increment_beyond_common[
        allocation_index[-1L]
      ]
    comparison[replicate, ] <-
      decomposition$comparisons$estimate[comparison_index]
    common_time_R2[replicate] <-
      decomposition$allocation$common_time_R2[1L]
    heterogeneity_increment_R2[replicate] <-
      decomposition$allocation$heterogeneity_increment_R2[1L]
    full_model_R2[replicate] <-
      decomposition$allocation$full_model_R2[1L]
  }
  list(
    allocated_R2 = allocated_R2,
    share_full = share_full,
    share_increment = share_increment,
    comparisons = comparison,
    common_time_R2 = common_time_R2,
    heterogeneity_increment_R2 = heterogeneity_increment_R2,
    full_model_R2 = full_model_R2,
    replicates = as.integer(replicates),
    seed = seed
  )
}

h02_dominance_quantile_interval <- function(x) {
  finite <- x[is.finite(x)]
  if (length(finite) == 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  stats::quantile(
    finite,
    probs = c(0.025, 0.975),
    names = FALSE
  ) |>
    stats::setNames(c("lower", "upper"))
}

h02_dominance_matrix_intervals <- function(x) {
  t(vapply(
    seq_len(ncol(x)),
    function(index) h02_dominance_quantile_interval(x[, index]),
    numeric(2)
  )) |>
    `rownames<-`(colnames(x))
}

h02_dominance_interval_method <- function(bootstrap) {
  paste(
    "95% percentile hierarchical cluster bootstrap of fixed predictions",
    "from all eight fitted subset models; sites, participants within sites,",
    "and participant-days within participants; conditional on the fitted",
    "subset models;",
    bootstrap$replicates,
    "replicates"
  )
}

h02_dominance_interval_summary <- function(
  point,
  bootstrap,
  run_id
) {
  component_names <- colnames(bootstrap$allocated_R2)
  allocation_interval <- h02_dominance_matrix_intervals(
    bootstrap$allocated_R2
  )
  share_full_interval <- h02_dominance_matrix_intervals(
    bootstrap$share_full
  )
  share_increment_interval <- h02_dominance_matrix_intervals(
    bootstrap$share_increment
  )
  common_interval <- h02_dominance_quantile_interval(
    bootstrap$common_time_R2
  )
  heterogeneity_interval <- h02_dominance_quantile_interval(
    bootstrap$heterogeneity_increment_R2
  )
  full_interval <- h02_dominance_quantile_interval(
    bootstrap$full_model_R2
  )
  component_index <- match(component_names, point$component)
  increment_lower <- rep(NA_real_, length(component_names))
  increment_upper <- rep(NA_real_, length(component_names))
  heterogeneity_components <- component_names != "common_time"
  increment_lower[heterogeneity_components] <-
    share_increment_interval[
      component_names[heterogeneity_components],
      "lower"
    ]
  increment_upper[heterogeneity_components] <-
    share_increment_interval[
      component_names[heterogeneity_components],
      "upper"
    ]
  tibble::tibble(
    run_id = run_id,
    component = component_names,
    allocated_R2 = point$allocated_R2[component_index],
    allocated_R2_lower_95 = allocation_interval[component_names, "lower"],
    allocated_R2_upper_95 = allocation_interval[component_names, "upper"],
    share_of_full_model_R2 = point$share_of_full_model_R2[component_index],
    share_of_full_model_R2_lower_95 = share_full_interval[
      component_names,
      "lower"
    ],
    share_of_full_model_R2_upper_95 = share_full_interval[
      component_names,
      "upper"
    ],
    share_of_increment_beyond_common = point$share_of_increment_beyond_common[
      component_index
    ],
    share_of_increment_beyond_common_lower_95 = increment_lower,
    share_of_increment_beyond_common_upper_95 = increment_upper,
    common_time_R2 = point$common_time_R2[component_index],
    common_time_R2_lower_95 = common_interval["lower"],
    common_time_R2_upper_95 = common_interval["upper"],
    heterogeneity_increment_R2 = point$heterogeneity_increment_R2[
      component_index
    ],
    heterogeneity_increment_R2_lower_95 = heterogeneity_interval["lower"],
    heterogeneity_increment_R2_upper_95 = heterogeneity_interval["upper"],
    full_model_R2 = point$full_model_R2[component_index],
    full_model_R2_lower_95 = full_interval["lower"],
    full_model_R2_upper_95 = full_interval["upper"],
    allocation_method = point$allocation_method[component_index],
    general_dominance_rank = point$general_dominance_rank[component_index],
    finite_bootstrap_replicates_allocated_R2 = colSums(is.finite(
      bootstrap$allocated_R2
    ))[component_names],
    finite_bootstrap_replicates_share_full = colSums(is.finite(
      bootstrap$share_full
    ))[component_names],
    finite_bootstrap_replicates_share_increment = c(
      common_time = NA_integer_,
      colSums(is.finite(bootstrap$share_increment))
    )[component_names],
    interval_method = h02_dominance_interval_method(bootstrap),
    bootstrap_seed = bootstrap$seed,
    bootstrap_replicates = bootstrap$replicates
  )
}

h02_dominance_comparison_interval_summary <- function(
  point,
  bootstrap,
  run_id
) {
  comparison_names <- colnames(bootstrap$comparisons)
  intervals <- h02_dominance_matrix_intervals(
    bootstrap$comparisons
  )
  point_index <- match(comparison_names, point$comparison_id)
  tibble::tibble(
    run_id = run_id,
    comparison_id = comparison_names,
    definition = point$definition[point_index],
    estimate = point$estimate[point_index],
    lower_95 = intervals[comparison_names, "lower"],
    upper_95 = intervals[comparison_names, "upper"],
    unit = point$unit[point_index],
    finite_bootstrap_replicates = colSums(is.finite(bootstrap$comparisons))[
      comparison_names
    ],
    interval_method = h02_dominance_interval_method(bootstrap),
    bootstrap_seed = bootstrap$seed,
    bootstrap_replicates = bootstrap$replicates
  )
}
