# Build H02-only model frames from prepared measurements.

h02_key <- c("site", "Id", "position", "local_date", "clock_bin")
h02_pair_key <- c("site", "Id", "local_date", "clock_bin")

h02_temporal_links <- function(root) {
  links <- readRDS(file.path(
    root,
    "results/intermediate/model_data/temporal_provenance/wall_outcome_links.rds"
  ))
  links |>
    dplyr::filter(.data$resolution == "30_minute") |>
    dplyr::transmute(
      site = as.character(.data$site),
      Id = as.character(.data$Id),
      position = as.character(.data$position),
      local_date = as.Date(.data$local_date),
      clock_bin = as.integer(.data$wall_bin_start_minute),
      source_bin_links = as.integer(.data$source_bin_links),
      source_utc_start = .data$source_utc_start,
      source_utc_end = .data$source_utc_end,
      relationship_type = as.character(.data$relationship_type),
      fold_provenance = as.character(.data$fold_provenance),
      spring_gap_provenance = as.character(.data$spring_gap_provenance),
      one_to_one_elapsed_coordinate = as.logical(
        .data$one_to_one_elapsed_coordinate
      )
    )
}

h02_assert_temporal_links <- function(links) {
  assert_unique_key(links, h02_key, "H02 30-minute temporal links")
  invalid_bin <- !links$clock_bin %in% seq.int(0L, 1410L, by = 30L)
  if (any(invalid_bin)) {
    h02_abort("Temporal links contain invalid 30-minute wall bins")
  }
  invisible(links)
}

h02_main_grid <- function(path, scenario_id, links) {
  grid <- readRDS(path)
  assert_columns(
    grid,
    c(
      h02_key,
      "metric_value_lx",
      "bin_admissible",
      "valid_medi_wall_minutes",
      "expected_wall_minutes",
      "failure_reason"
    ),
    "main H02 grid"
  )
  assert_unique_key(grid, h02_key, "main H02 grid")
  if (
    any(grid$expected_wall_minutes != 30L) ||
      any(!grid$clock_bin %in% seq.int(0L, 1410L, by = 30L))
  ) {
    h02_abort("Main H02 input does not have exact 30-minute wall-clock support")
  }
  support_ok <- grid$bin_admissible ==
    (grid$valid_medi_wall_minutes >= 15L &
      is.finite(grid$metric_value_lx))
  if (any(!support_ok)) {
    h02_abort(
      "Main H02 input violates the >=15 valid-minute support rule"
    )
  }
  dplyr::transmute(
    grid,
    data_scenario_id = scenario_id,
    site = as.character(.data$site),
    Id = as.character(.data$Id),
    position = as.character(.data$position),
    local_date = as.Date(.data$local_date),
    clock_bin = as.integer(.data$clock_bin),
    metric_value_lx = as.numeric(.data$metric_value_lx),
    bin_admissible = as.logical(.data$bin_admissible),
    valid_medi_wall_minutes = as.integer(.data$valid_medi_wall_minutes),
    support_available = TRUE,
    failure_reason = as.character(.data$failure_reason)
  ) |>
    left_join_checked(
      links,
      by = h02_key,
      relationship = "one-to-one",
      x_name = "main H02 grid",
      y_name = "temporal links"
    )
}

h02_clock_bin_from_proxy <- function(x) {
  lt <- as.POSIXlt(x, tz = "UTC")
  as.integer(lt$hour * 60L + lt$min)
}

h02_alternative_grid <- function(path, links) {
  grid <- readRDS(path)
  assert_columns(
    grid,
    c(
      "scenario_id",
      "position",
      "site",
      "Id",
      "local_date",
      "local_clock_datetime_utc_proxy",
      "medi_arithmetic_mean_lx"
    ),
    "alternative-preprocessing H02 grid"
  )
  grid <- dplyr::transmute(
    grid,
    data_scenario_id = "alternative_preprocessing",
    site = as.character(.data$site),
    Id = as.character(.data$Id),
    position = as.character(.data$position),
    local_date = as.Date(.data$local_date),
    clock_bin = h02_clock_bin_from_proxy(
      .data$local_clock_datetime_utc_proxy
    ),
    metric_value_lx = as.numeric(.data$medi_arithmetic_mean_lx),
    bin_admissible = is.finite(.data$medi_arithmetic_mean_lx),
    valid_medi_wall_minutes = NA_integer_,
    support_available = FALSE,
    failure_reason = dplyr::if_else(
      is.finite(.data$medi_arithmetic_mean_lx),
      NA_character_,
      "alternative_preprocessing_nonfinite"
    )
  )
  assert_unique_key(grid, h02_key, "alternative-preprocessing H02 grid")
  if (any(!grid$clock_bin %in% seq.int(0L, 1410L, by = 30L))) {
    h02_abort("Alternative-preprocessing H02 input contains non-30-minute bins")
  }
  left_join_checked(
    grid,
    links,
    by = h02_key,
    relationship = "one-to-one",
    x_name = "alternative-preprocessing H02 grid",
    y_name = "temporal links"
  )
}

h02_assign_sequences <- function(frame) {
  frame <- frame |>
    dplyr::arrange(
      .data$site,
      .data$Id,
      .data$local_date,
      .data$source_utc_start,
      .data$clock_bin
    ) |>
    dplyr::group_by(.data$site, .data$Id, .data$local_date) |>
    dplyr::mutate(
      previous_source_utc_end = dplyr::lag(.data$source_utc_end),
      previous_one_to_one = dplyr::lag(
        .data$one_to_one_elapsed_coordinate
      ),
      elapsed_from_previous_seconds = as.numeric(difftime(
        .data$source_utc_start,
        .data$previous_source_utc_end,
        units = "secs"
      )),
      ar_start_reason = dplyr::case_when(
        dplyr::row_number() == 1L ~ "participant_day_start",
        !.data$one_to_one_elapsed_coordinate ~ "non_one_to_one_wall_outcome",
        !dplyr::coalesce(.data$previous_one_to_one, TRUE) ~
          "after_non_one_to_one_wall_outcome",
        is.na(.data$elapsed_from_previous_seconds) ~
          "missing_elapsed_coordinate",
        .data$elapsed_from_previous_seconds != 0 ~ "elapsed_discontinuity",
        TRUE ~ "continuous"
      ),
      AR_start = .data$ar_start_reason != "continuous",
      sequence_number = cumsum(.data$AR_start)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      participant_key = paste(.data$site, .data$Id, sep = "::"),
      participant_day_key = paste(
        .data$site,
        .data$Id,
        .data$local_date,
        sep = "::"
      ),
      true_elapsed_sequence_id = paste(
        .data$participant_day_key,
        sprintf("S%03d", .data$sequence_number),
        sep = "::"
      ),
      time_hour = (.data$clock_bin + 15) / 60,
      response = h02_transform(.data$metric_value_lx)
    ) |>
    dplyr::select(
      -"previous_source_utc_end",
      -"previous_one_to_one",
      -"sequence_number"
    )

  if (
    anyNA(frame$source_utc_start) ||
      anyNA(frame$source_utc_end) ||
      anyNA(frame$one_to_one_elapsed_coordinate)
  ) {
    h02_abort(
      "A finite H02 outcome lacks its true-elapsed-time mapping"
    )
  }
  if (
    any(!frame$clock_bin %in% seq.int(0L, 1410L, by = 30L)) ||
      any(frame$time_hour <= 0 | frame$time_hour >= 24)
  ) {
    h02_abort("H02 fitted frame violates 30-minute midpoint support")
  }
  first_by_day <- !duplicated(frame[c("site", "Id", "local_date")])
  if (!all(frame$AR_start[first_by_day])) {
    h02_abort("AR(1) does not reset at every participant-day")
  }
  if (any(frame$elapsed_from_previous_seconds[!frame$AR_start] != 0)) {
    h02_abort("AR(1) crosses a true-elapsed discontinuity")
  }
  if (
    any(
      !frame$one_to_one_elapsed_coordinate[!frame$AR_start] |
        !dplyr::lag(
          frame$one_to_one_elapsed_coordinate,
          default = TRUE
        )[!frame$AR_start]
    )
  ) {
    h02_abort("AR(1) crosses a non-one-to-one fall-back outcome")
  }
  frame
}

h02_model_frame <- function(grid) {
  fitted <- grid |>
    dplyr::filter(
      .data$bin_admissible,
      is.finite(.data$metric_value_lx)
    )
  if (nrow(fitted) == 0L) {
    h02_abort("H02 grid has no admissible finite outcomes")
  }
  h02_assign_sequences(fitted)
}

h02_paired_frames <- function(glasses, chest) {
  glasses_keys <- dplyr::distinct(
    glasses,
    dplyr::across(
      dplyr::all_of(h02_pair_key)
    )
  )
  chest_keys <- dplyr::distinct(
    chest,
    dplyr::across(
      dplyr::all_of(h02_pair_key)
    )
  )
  common <- dplyr::inner_join(
    glasses_keys,
    chest_keys,
    by = h02_pair_key,
    relationship = "one-to-one"
  )
  if (nrow(common) == 0L) {
    h02_abort("No common near-eye/chest H02 bins")
  }
  list(
    glasses = glasses |>
      dplyr::semi_join(common, by = h02_pair_key) |>
      h02_assign_sequences(),
    chest = chest |>
      dplyr::semi_join(common, by = h02_pair_key) |>
      h02_assign_sequences()
  )
}

h02_sample_counts <- function(frame, run_id) {
  overall <- tibble::tibble(
    run_id = run_id,
    site = "ALL_SITES",
    participants = dplyr::n_distinct(frame$participant_key),
    participant_days = dplyr::n_distinct(frame$participant_day_key),
    observations_30_minute = nrow(frame),
    sites = dplyr::n_distinct(frame$site),
    exact_zero_observations = sum(frame$metric_value_lx == 0),
    ar_sequences = dplyr::n_distinct(frame$true_elapsed_sequence_id)
  )
  by_site <- frame |>
    dplyr::group_by(.data$site) |>
    dplyr::summarise(
      participants = dplyr::n_distinct(.data$participant_key),
      participant_days = dplyr::n_distinct(.data$participant_day_key),
      observations_30_minute = dplyr::n(),
      sites = 1L,
      exact_zero_observations = sum(.data$metric_value_lx == 0),
      ar_sequences = dplyr::n_distinct(.data$true_elapsed_sequence_id),
      .groups = "drop"
    ) |>
    dplyr::mutate(run_id = run_id, .before = 1L)
  dplyr::bind_rows(overall, by_site)
}

h02_support_audit <- function(grid, scenario_id, placement) {
  tibble::tibble(
    data_scenario_id = scenario_id,
    placement = placement,
    grid_rows = nrow(grid),
    participants = dplyr::n_distinct(grid$Id),
    participant_days = dplyr::n_distinct(
      paste(grid$site, grid$Id, grid$local_date)
    ),
    sites = dplyr::n_distinct(grid$site),
    finite_admissible = sum(
      grid$bin_admissible & is.finite(grid$metric_value_lx)
    ),
    unsupported_or_nonfinite = sum(
      !grid$bin_admissible | !is.finite(grid$metric_value_lx)
    ),
    support_available = all(grid$support_available),
    minimum_valid_minutes_admissible = if (all(grid$support_available)) {
      min(grid$valid_medi_wall_minutes[grid$bin_admissible])
    } else {
      NA_integer_
    },
    maximum_valid_minutes_inadmissible = if (all(grid$support_available)) {
      max(grid$valid_medi_wall_minutes[!grid$bin_admissible])
    } else {
      NA_integer_
    },
    non_one_to_one_finite = sum(
      grid$bin_admissible &
        is.finite(grid$metric_value_lx) &
        !grid$one_to_one_elapsed_coordinate
    )
  )
}
