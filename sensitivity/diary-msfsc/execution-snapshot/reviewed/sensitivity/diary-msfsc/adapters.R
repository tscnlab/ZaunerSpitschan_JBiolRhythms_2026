# Generalize report calculations to changed samples and any number of windows.
# Model fitting and covariance calculations remain in the frozen source snapshot.

diary_daylight_design <- function(d, specification, prefix) {
  ref <- reference_days(d)
  probs <- switch(
    specification,
    fixed = rep(NA_real_, 3),
    quartiles = c(.25, .5, .75),
    p10_p90 = c(.1, .5, .9)
  )
  targets <- tibble(
    scenario = c("Shorter", "Typical", "Longer"),
    probability = probs,
    photo = if (specification == "fixed") c(-1, 0, 1) else
      daylight_weighted_quantile(ref$photo, ref$w, probs)
  )
  support <- bind_rows(lapply(seq_len(nrow(targets)), function(i) {
    d |>
      filter(abs(photo - targets$photo[i]) <= 1, quartile %in% c("Q1", "Q4")) |>
      mutate(
        quartile = as.character(quartile),
        daytype = as.character(daytype)
      ) |>
      summarise(
        participants = n_distinct(participant_key),
        participant_days = n_distinct(participant_day_key),
        half_hours = n(),
        .by = c(quartile, daytype)
      ) |>
      complete(
        quartile = c("Q1", "Q4"),
        daytype = c("Work", "Free"),
        fill = list(participants = 0L, participant_days = 0L, half_hours = 0L)
      ) |>
      mutate(
        scenario = targets$scenario[i],
        photo_target = targets$photo[i],
        photo_tolerance_hours = 1
      )
  }))
  site <- d |>
    filter(quartile %in% c("Q1", "Q4")) |>
    mutate(quartile = as.character(quartile)) |>
    summarise(
      participants = n_distinct(participant_key),
      participant_days = n_distinct(participant_day_key),
      half_hours = n(),
      .by = c(site, quartile)
    ) |>
    complete(
      site,
      quartile = c("Q1", "Q4"),
      fill = list(participants = 0L, participant_days = 0L, half_hours = 0L)
    )
  design <- list(
    targets = targets,
    free_weight = sum(ref$w[ref$daytype == "Free"]),
    support = support,
    site = site,
    formula = deparse(daylight_formula()),
    version = daylight_version,
    reporting_version = daylight_reporting_version,
    specification = specification,
    prefix = prefix,
    support_n = daylight_support_n
  )
  save_table(targets, paste0(prefix, "_scenario_targets"))
  save_table(support, paste0(prefix, "_scenario_samples"))
  save_table(site, paste0(prefix, "_site_composition"))
  atomic_rds(design, paste0("results/data/", prefix, "_design.rds"))
  design
}

diary_daylight_intervals <- function(bundle, inference, design, windows) {
  mm <- daylight_matrices(bundle, design)
  covariance <- covariance_info(bundle$fit)
  rows <- list()
  for (i in seq_len(nrow(windows)))
    for (id in names(mm$D)) {
      ix <- which(
        wrap24(mm$x - windows$start[i]) < windows$duration_hours[i] - 1e-10
      )
      z <- linear_summary(
        matrix(colMeans(mm$D[[id]][ix, , drop = FALSE]), nrow = 1),
        coef(bundle$fit),
        covariance$V
      )
      support <- filter(inference$contrasts, contrast == id, x %in% mm$x[ix])
      rows[[length(rows) + 1L]] <- bind_cols(
        windows[i, c("window_source", "start", "end", "duration_hours")],
        tibble(
          contrast = id,
          phase_label = ifelse(windows$end[i] < 0, "Pre-MSFsc", "Post-MSFsc")
        ),
        z
      ) |>
        mutate(
          ratio = 10^estimate,
          ratio_lower = 10^lower,
          ratio_upper = 10^upper,
          supported = all(support$supported),
          minimum_cell_participants = min(support$minimum_cell_participants),
          grid_points = length(ix),
          covariance = covariance$label,
          interval_method = "Mean log contrast; fixed window; unadjusted 95% Gaussian interval"
        )
    }
  out <- bind_rows(rows)
  save_table(out, paste0(design$prefix, "_interval_summaries"))
  out
}

diary_solar_summary <- function(bundle, windows) {
  d <- bundle$data
  context <- readRDS(input_paths()[["glasses_context"]])
  keys <- c("site", "Id", "position", "local_date", "clock_bin")
  solar <- select(
    context,
    all_of(keys),
    civil_dawn_wall_minute,
    civil_dusk_wall_minute,
    solar_depression_deg
  )
  stopifnot(!anyDuplicated(solar[keys]))
  joined <- left_join(d, solar, by = keys, relationship = "one-to-one")
  stopifnot(
    nrow(joined) == nrow(d),
    !anyNA(joined$civil_dawn_wall_minute),
    !anyNA(joined$civil_dusk_wall_minute),
    all(joined$solar_depression_deg == 6)
  )
  days <- joined |>
    distinct(
      participant_key,
      participant_day_key,
      site,
      quartile,
      msf_hours,
      photoperiod_hours,
      civil_dawn_wall_minute,
      civil_dusk_wall_minute
    )
  ref <- days |> group_by(participant_key) |> mutate(w = 1 / n()) |> ungroup()
  split <- daylight_weighted_quantile(ref$photoperiod_hours, ref$w, .5)
  days <- mutate(
    days,
    daylength = if_else(photoperiod_hours <= split, "Shorter", "Longer")
  )
  extreme <- filter(days, quartile %in% c("Q1", "Q4")) |>
    group_by(quartile, daylength, participant_key) |>
    mutate(w = 1 / n()) |>
    ungroup() |>
    group_by(quartile, daylength) |>
    mutate(w = w / n_distinct(participant_key)) |>
    ungroup()
  x <- seq(-12, 12 - cfg$step, by = cfg$step)
  estimates <- bind_rows(lapply(x, function(xx) {
    extreme |>
      mutate(
        clock_minute = wrap24(xx + msf_hours) * 60,
        available = ((clock_minute - civil_dawn_wall_minute) %% 1440) <=
          ((civil_dusk_wall_minute - civil_dawn_wall_minute) %% 1440)
      ) |>
      summarise(
        availability = sum(w * available),
        participants = n_distinct(participant_key),
        participant_days = n_distinct(participant_day_key),
        .by = c(quartile, daylength)
      ) |>
      mutate(x = xx)
  }))
  interval_summaries <- bind_rows(lapply(seq_len(nrow(windows)), function(i) {
    estimates |>
      filter(
        wrap24(x - windows$start[i]) < windows$duration_hours[i] - 1e-10
      ) |>
      summarise(
        availability_percent = 100 * mean(availability),
        grid_points = n(),
        participants = first(participants),
        participant_days = first(participant_days),
        .by = c(quartile, daylength)
      ) |>
      mutate(
        window_source = windows$window_source[i],
        start = windows$start[i],
        end = windows$end[i]
      )
  }))
  save_table(estimates, "daylight_civil_availability")
  save_table(interval_summaries, "daylight_civil_interval_summaries")
  save_table(
    tibble(weighted_median_hours = split, solar_depression_degrees = 6),
    "daylight_civil_definition"
  )
  invisible(interval_summaries)
}
