# Compare diary midsleep with questionnaire MSFsc, preserving waking-day labels.

## ---- sleep-diary-inputs
source(".Rprofile")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
sleep_diary_manifest <- read.csv("input-data/manifest.csv") |>
  filter(role %in% c("sleepdiaries", "chronotype"))
stopifnot(identical(
  unname(tools::md5sum(sleep_diary_manifest$path)),
  sleep_diary_manifest$md5
))
sleep_diary_raw <- readRDS(
  "input-data/model_data/normalized_inputs/sleepdiaries.rds"
)
sleep_diary_chronotype <- readRDS(
  "input-data/model_data/normalized_inputs/chronotype.rds"
) |>
  transmute(site, Id, msfsc_hours = as.numeric(msf_sc) / 3600)
sleep_diary_primary <- readRDS("results/data/prepared.rds")$frames$glasses |>
  distinct(site, Id, local_date)
stopifnot(!anyDuplicated(sleep_diary_chronotype[c("site", "Id")]))

## ---- sleep-diary-prepare
sleep_diary_records <- sleep_diary_raw |>
  left_join(
    sleep_diary_chronotype,
    by = c("site", "Id"),
    relationship = "many-to-one"
  ) |>
  left_join(
    sleep_diary_primary |>
      distinct(site, Id) |>
      mutate(in_primary_sample = TRUE),
    by = c("site", "Id"),
    relationship = "many-to-one"
  ) |>
  mutate(
    participant = paste(site, Id, sep = "::"),
    in_primary_sample = coalesce(in_primary_sample, FALSE),
    waking_date = as.Date(wake_wall, tz = "UTC"),
    daytype = factor(
      as.character(daytype2),
      levels = c("a work day", "a free day"),
      labels = c("Work", "Free")
    ),
    sleep_period_hours = as.numeric(difftime(
      wake_wall,
      sleep_onset_wall,
      units = "hours"
    )),
    midpoint_wall = sleep_onset_wall + sleep_period_hours * 3600 / 2,
    midsleep_hours = (as.numeric(midpoint_wall) / 3600) %% 24,
    # A fixed noon cut keeps times either side of midnight adjacent.
    midsleep_unwrapped_hours = ((midsleep_hours + 12) %% 24) - 12,
    msfsc_unwrapped_hours = ((msfsc_hours + 12) %% 24) - 12,
    difference_minutes = (((midsleep_hours - msfsc_hours + 12) %% 24) - 12) *
      60,
    dst_crossing = sleep_onset_utc_offset_minutes != wake_utc_offset_minutes,
    exclusion = case_when(
      !is.finite(msfsc_hours) ~ "Missing MSFsc",
      is.na(sleep_onset_wall) | is.na(wake_wall) ~
        "Missing sleep onset or awakening",
      !sleep_interval_analysis_eligible | sleep_interval_quarantined ~
        "Upstream interval quarantine",
      sleep_period_hours <= 0 | sleep_period_hours >= 24 ~
        "Invalid sleep-period duration",
      is.na(daytype) ~ "Missing waking-day type",
      TRUE ~ "Included"
    )
  ) |>
  left_join(
    sleep_diary_primary |>
      rename(waking_date = local_date) |>
      mutate(waking_date_in_light_analysis = TRUE),
    by = c("site", "Id", "waking_date"),
    relationship = "many-to-one"
  ) |>
  mutate(
    waking_date_in_light_analysis = coalesce(
      waking_date_in_light_analysis,
      FALSE
    )
  )
sleep_diary_valid <- sleep_diary_records |> filter(exclusion == "Included")
sleep_diary_nights <- sleep_diary_valid |> filter(in_primary_sample)
sleep_diary_flow <- bind_rows(
  sleep_diary_records |>
    summarise(
      stage = "All imported diary entries",
      participants = n_distinct(participant),
      nights = n()
    ),
  sleep_diary_records |>
    filter(in_primary_sample) |>
    summarise(
      stage = "Entries from primary exposure participants",
      participants = n_distinct(participant),
      nights = n()
    ),
  sleep_diary_nights |>
    summarise(
      stage = "Valid primary-sample diary nights",
      participants = n_distinct(participant),
      nights = n()
    ),
  sleep_diary_nights |>
    filter(waking_date_in_light_analysis) |>
    summarise(
      stage = "Also have waking-day light observations",
      participants = n_distinct(participant),
      nights = n()
    ),
  sleep_diary_valid |>
    summarise(
      stage = "Valid nights in all available participants",
      participants = n_distinct(participant),
      nights = n()
    )
)
sleep_diary_flow

## ---- sleep-diary-correlation-helper
# Resample whole people, carrying all their nights together.
sleep_diary_correlations <- function(nights, context, scope, draws = 2000L) {
  selected <- if (context == "Overall") nights else
    filter(nights, daytype == context)
  people <- selected |>
    summarise(
      msfsc = first(msfsc_unwrapped_hours),
      midsleep = mean(midsleep_unwrapped_hours),
      .by = participant
    )
  stopifnot(nrow(people) >= 3L, stats::sd(people$msfsc) > 0)
  night_indices <- split(seq_len(nrow(selected)), selected$participant)
  boot <- replicate(draws, {
    chosen <- sample.int(nrow(people), nrow(people), replace = TRUE)
    ids <- people$participant[chosen]
    rows <- unlist(night_indices[ids], use.names = FALSE)
    c(
      night = stats::cor(
        selected$msfsc_unwrapped_hours[rows],
        selected$midsleep_unwrapped_hours[rows]
      ),
      participant_mean = stats::cor(
        people$msfsc[chosen],
        people$midsleep[chosen]
      )
    )
  })
  stopifnot(all(is.finite(boot)))
  tibble::tibble(
    scope = scope,
    context = context,
    participants = nrow(people),
    nights = nrow(selected),
    night_r = stats::cor(
      selected$msfsc_unwrapped_hours,
      selected$midsleep_unwrapped_hours
    ),
    night_lower = unname(stats::quantile(boot["night", ], 0.025)),
    night_upper = unname(stats::quantile(boot["night", ], 0.975)),
    participant_mean_r = stats::cor(people$msfsc, people$midsleep),
    participant_mean_lower = unname(stats::quantile(
      boot["participant_mean", ],
      0.025
    )),
    participant_mean_upper = unname(stats::quantile(
      boot["participant_mean", ],
      0.975
    )),
    median_difference_minutes = stats::median(selected$difference_minutes),
    median_absolute_difference_minutes = stats::median(abs(
      selected$difference_minutes
    )),
    bootstrap_draws = draws
  )
}

## ---- sleep-diary-correlations
sleep_diary_correlation_table <- withr::with_seed(20261106L, {
  rows <- list()
  for (scope in c(
    "Primary exposure participants",
    "All available participants"
  )) {
    nights <- if (scope == "Primary exposure participants")
      sleep_diary_nights else sleep_diary_valid
    for (context in c("Overall", "Work", "Free")) {
      rows[[length(rows) + 1L]] <- sleep_diary_correlations(
        nights,
        context,
        scope
      )
    }
  }
  bind_rows(rows)
})
sleep_diary_correlation_table |>
  filter(scope == "Primary exposure participants") |>
  select(
    context,
    participants,
    nights,
    night_r,
    night_lower,
    night_upper,
    participant_mean_r
  )

## ---- sleep-diary-variability
sleep_diary_contexts <- bind_rows(
  sleep_diary_nights |> mutate(context = "Overall"),
  sleep_diary_nights |> mutate(context = as.character(daytype))
)
sleep_diary_participant_sd <- sleep_diary_contexts |>
  summarise(
    nights = n(),
    msfsc_hours = first(msfsc_hours),
    mean_midsleep_hours = mean(midsleep_unwrapped_hours),
    sd_minutes = if (n() > 2L) stats::sd(midsleep_unwrapped_hours) * 60 else
      NA_real_,
    .by = c(site, Id, participant, context)
  ) |>
  mutate(eligible_for_sd = nights > 2L)
sleep_diary_sd_summary <- sleep_diary_participant_sd |>
  filter(eligible_for_sd) |>
  summarise(
    participants = n(),
    median_sd_minutes = stats::median(sd_minutes),
    q1_sd_minutes = unname(stats::quantile(sd_minutes, 0.25)),
    q3_sd_minutes = unname(stats::quantile(sd_minutes, 0.75)),
    minimum_nights = min(nights),
    maximum_nights = max(nights),
    nights = sum(nights),
    .by = context
  ) |>
  arrange(match(context, c("Overall", "Work", "Free")))
sleep_diary_sd_summary

## ---- sleep-diary-msfsc
sleep_diary_daytype_means <- sleep_diary_nights |>
  summarise(
    nights = n(),
    midsleep = mean(midsleep_unwrapped_hours),
    sleep_duration = mean(sleep_period_hours),
    .by = c(site, Id, participant, msfsc_unwrapped_hours, daytype)
  )
sleep_diary_msfsc_participants <- sleep_diary_daytype_means |>
  pivot_wider(
    names_from = daytype,
    values_from = c(nights, midsleep, sleep_duration),
    values_fill = list(nights = 0L)
  ) |>
  rename(
    work_nights = nights_Work,
    free_nights = nights_Free,
    msw_hours = midsleep_Work,
    msf_hours = midsleep_Free,
    sd_work_hours = sleep_duration_Work,
    sd_free_hours = sleep_duration_Free
  ) |>
  mutate(
    eligible = work_nights > 0L & free_nights > 0L,
    mean_sleep_duration_hours = if_else(
      eligible,
      (work_nights * sd_work_hours + free_nights * sd_free_hours) /
        (work_nights + free_nights),
      NA_real_
    ),
    correction_hours = if_else(
      eligible,
      pmax(0, sd_free_hours - mean_sleep_duration_hours) / 2,
      NA_real_
    ),
    diary_msfsc_hours = msf_hours - correction_hours,
    difference_minutes = (diary_msfsc_hours - msfsc_unwrapped_hours) * 60
  )
sleep_diary_msfsc_eligible <- sleep_diary_msfsc_participants |>
  filter(eligible)
sleep_diary_msfsc_flow <- sleep_diary_msfsc_participants |>
  summarise(
    primary_participants = n(),
    eligible_participants = sum(eligible),
    missing_work_nights = sum(work_nights == 0L),
    missing_free_nights = sum(free_nights == 0L),
    work_nights = sum(work_nights[eligible]),
    free_nights = sum(free_nights[eligible]),
    corrected_participants = sum(correction_hours > 0, na.rm = TRUE),
    median_correction_minutes = stats::median(correction_hours, na.rm = TRUE) *
      60
  )
sleep_diary_msfsc_flow

## ---- sleep-diary-msfsc-summary
sleep_diary_msfsc_summary <- withr::with_seed(20261107L, {
  people <- sleep_diary_msfsc_eligible
  metrics <- c("msw_hours", "msf_hours", "diary_msfsc_hours")
  draws <- 2000L
  boot <- replicate(draws, {
    rows <- sample.int(nrow(people), nrow(people), replace = TRUE)
    vapply(
      metrics,
      function(metric) {
        stats::cor(people$msfsc_unwrapped_hours[rows], people[[metric]][rows])
      },
      numeric(1)
    )
  })
  stopifnot(all(is.finite(boot)))
  bind_rows(lapply(metrics, function(metric) {
    difference <- (people[[metric]] - people$msfsc_unwrapped_hours) * 60
    tibble::tibble(
      metric = metric,
      participants = nrow(people),
      r = stats::cor(people$msfsc_unwrapped_hours, people[[metric]]),
      lower = unname(stats::quantile(boot[metric, ], 0.025)),
      upper = unname(stats::quantile(boot[metric, ], 0.975)),
      median_difference_minutes = stats::median(difference),
      q1_difference_minutes = unname(stats::quantile(difference, 0.25)),
      q3_difference_minutes = unname(stats::quantile(difference, 0.75)),
      median_absolute_difference_minutes = stats::median(abs(difference)),
      bootstrap_draws = draws
    )
  }))
})
sleep_diary_msfsc_summary

## ---- sleep-diary-scatter-data
sleep_diary_scatter_nights_data <- sleep_diary_nights |>
  select(
    site,
    Id,
    participant,
    waking_date,
    daytype,
    msfsc_hours,
    midsleep_hours,
    msfsc_unwrapped_hours,
    midsleep_unwrapped_hours,
    sleep_period_hours,
    difference_minutes,
    dst_crossing,
    waking_date_in_light_analysis
  )
sleep_diary_scatter_means_data <- sleep_diary_scatter_nights_data |>
  summarise(
    nights = n(),
    msfsc_unwrapped_hours = first(msfsc_unwrapped_hours),
    midsleep_unwrapped_hours = mean(midsleep_unwrapped_hours),
    .by = c(site, Id, participant, daytype)
  )
sleep_diary_scatter_data <- bind_rows(
  sleep_diary_scatter_nights_data |>
    mutate(panel = "A  Individual nights", nights = 1L),
  sleep_diary_scatter_means_data |>
    mutate(panel = "B  Participant means"),
  sleep_diary_msfsc_eligible |>
    transmute(
      site,
      Id,
      participant,
      msfsc_unwrapped_hours,
      midsleep_unwrapped_hours = diary_msfsc_hours,
      nights = work_nights + free_nights,
      panel = "C  Diary-derived MSFsc"
    )
)
sleep_diary_scatter_correlations_data <- bind_rows(
  sleep_diary_correlation_table |>
    mutate(panel = "A  Individual nights", r = night_r, n = nights),
  sleep_diary_correlation_table |>
    mutate(panel = "B  Participant means", r = participant_mean_r, n = participants)
) |>
  filter(scope == "Primary exposure participants") |>
  transmute(
    panel,
    daytype = factor(context, levels = c("Overall", "Work", "Free")),
    r,
    n,
    label = sprintf("%s: r = %.2f; n = %d", context, r, n),
    x = 0.2,
    y = 11.6 - 0.7 * (as.integer(daytype) - 1)
  )
sleep_diary_scatter_correlations_data <- bind_rows(
  sleep_diary_scatter_correlations_data,
  sleep_diary_msfsc_summary |>
    filter(metric == "diary_msfsc_hours") |>
    transmute(
      panel = "C  Diary-derived MSFsc",
      r,
      n = participants,
      label = sprintf("r = %.2f; n = %d", r, n),
      x = 0.2,
      y = 11.6
    )
)
# The marginal histogram counts participant means, with common 30-minute bins.
sleep_diary_scatter_histogram_data <- sleep_diary_scatter_means_data |>
  mutate(bin_lower = floor(midsleep_unwrapped_hours * 2) / 2) |>
  count(daytype, bin_lower, name = "participants") |>
  complete(
    daytype,
    bin_lower = seq(0, 11.5, by = 0.5),
    fill = list(participants = 0L)
  ) |>
  mutate(
    bin_upper = bin_lower + 0.5,
    bar_lower = bin_lower + if_else(daytype == "Free", 0.25, 0),
    bar_upper = bar_lower + 0.25
  )

## ---- sleep-diary-plot
# Solarized accents from ggthemr's palette, used directly without a theme dependency.
# https://github.com/Mikata-Project/ggthemr/blob/master/R/palettes.R
sleep_diary_day_colors <- c(Work = "#268bd2", Free = "#dc322f")
sleep_diary_scatter_correlations_data <- sleep_diary_scatter_correlations_data |>
  mutate(
    label_color = coalesce(
      sleep_diary_day_colors[as.character(daytype)],
      "black"
    )
  )
sleep_diary_panels <- list()
for (panel in unique(sleep_diary_scatter_data$panel)) {
  plot_data <- filter(sleep_diary_scatter_data, .data$panel == .env$panel)
  corrected <- panel == "C  Diary-derived MSFsc"
  p <- ggplot(plot_data, aes(msfsc_unwrapped_hours, midsleep_unwrapped_hours)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.55
    )
  if (corrected) {
    p <- p +
      geom_point(
        shape = 21,
        fill = "#2aa198B3",
        color = "black",
        size = 2.1,
        stroke = 0.3
      )
  } else {
    p <- p +
      geom_point(
        aes(fill = daytype, shape = daytype),
        color = "black",
        size = 2.1,
        stroke = 0.3
      ) +
      scale_fill_manual(values = scales::alpha(sleep_diary_day_colors, 0.7)) +
      scale_shape_manual(values = c(Work = 21, Free = 24)) +
      labs(fill = "Waking day", shape = "Waking day")
  }
  p <- p +
    geom_text(
      data = filter(
        sleep_diary_scatter_correlations_data,
        .data$panel == .env$panel
      ),
      aes(x = x, y = y, label = label, color = label_color),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.5
    ) +
    scale_color_identity() +
    scale_x_continuous(
      breaks = seq(0, 12, 3),
      labels = sprintf("%02d:00", seq(0, 12, 3))
    ) +
    scale_y_continuous(
      breaks = seq(0, 12, 3),
      labels = sprintf("%02d:00", seq(0, 12, 3))
    ) +
    coord_cartesian(xlim = c(0, 12), ylim = c(0, 12), expand = FALSE) +
    labs(
      title = sub("^[A-C] +", "", panel),
      tag = substr(panel, 1, 1),
      x = "Questionnaire MSFsc",
      y = if (corrected) "Diary MSFsc equivalent" else "Diary midsleep"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12),
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = "topleft",
      plot.margin = margin(
        6,
        if (panel == "B  Participant means") 0 else 12,
        6,
        6
      )
    )
  sleep_diary_panels[[panel]] <- p
}
sleep_diary_marginal <- ggplot(sleep_diary_scatter_histogram_data) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = participants,
      ymin = bar_lower,
      ymax = bar_upper,
      fill = daytype
    ),
    color = "black",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = scales::alpha(sleep_diary_day_colors, 0.7)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = seq(0, 12, 3)) +
  coord_cartesian(ylim = c(0, 12), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
sleep_diary_scatter <- patchwork::wrap_plots(
  sleep_diary_panels[[1]],
  sleep_diary_panels[[2]],
  sleep_diary_marginal,
  sleep_diary_panels[[3]],
  nrow = 1,
  widths = c(1, 1, 0.25, 1),
  guides = "collect"
) &
  theme(legend.position = "none")
if (isTRUE(getOption("knitr.in.progress"))) print(sleep_diary_scatter)

## ---- sleep-diary-checks
sleep_diary_dst_audit <- sleep_diary_nights |>
  mutate(
    midpoint_utc = sleep_onset_utc +
      as.numeric(difftime(wake_utc, sleep_onset_utc, units = "secs")) / 2
  )
sleep_diary_dst_audit$elapsed_midpoint_local_label <- vapply(
  seq_len(nrow(sleep_diary_dst_audit)),
  function(i)
    format(
      sleep_diary_dst_audit$midpoint_utc[i],
      tz = sleep_diary_dst_audit$site_timezone[i],
      format = "%Y-%m-%d %H:%M:%S"
    ),
  character(1)
)
sleep_diary_dst_audit <- sleep_diary_dst_audit |>
  mutate(
    elapsed_midpoint_wall = as.POSIXct(
      elapsed_midpoint_local_label,
      tz = "UTC"
    ),
    elapsed_midpoint_hours = (as.numeric(elapsed_midpoint_wall) / 3600) %% 24,
    elapsed_minus_wall_minutes = (((elapsed_midpoint_hours -
      midsleep_hours +
      12) %%
      24) -
      12) *
      60
  )
sleep_diary_dst_sensitivity <- sleep_diary_dst_audit |>
  summarise(
    night_r_wall = stats::cor(msfsc_unwrapped_hours, midsleep_unwrapped_hours),
    night_r_elapsed = stats::cor(
      msfsc_unwrapped_hours,
      ((elapsed_midpoint_hours + 12) %% 24) - 12
    ),
    nights_with_different_midpoint = sum(abs(elapsed_minus_wall_minutes) > 1e-6)
  )
sleep_diary_valid_day_keys <- sleep_diary_records |> filter(!is.na(waking_date))
sleep_diary_checks <- tibble::tibble(
  check = c(
    "One diary record per participant and waking date",
    "Primary participants match the exposure cohort",
    "Sleep onset equals attempted sleep plus latency",
    "All valid durations are positive and shorter than 24 hours",
    "Midpoints use dated timestamps across midnight",
    "Midnight-adjacent clocks remain adjacent on the noon branch",
    "No observed noon-branch discontinuity",
    "Diary day types are retained on their waking-day records",
    "SD requires at least three nights within its context",
    "SD is unavailable at two or fewer nights",
    "Elapsed and wall midpoint agree away from DST changes",
    "All bootstrap confidence limits are finite",
    "SD summaries retain per-participant night counts",
    "Panel B has one mean per participant and day type",
    "Panel B includes every plotted night in its day-type mean",
    "Panel B correlations match the participant-mean results",
    "Diary correction requires both work and free nights",
    "Count-weighted duration equals the average of the retained nights",
    "Correction is zero when free-day duration does not exceed work-day duration",
    "Correction is nonnegative and advances or retains free-day midsleep",
    "Corrected correlation confidence limits are finite",
    "Panel C has one corrected value per eligible participant",
    "Marginal histogram counts every participant-day-type mean",
    "All plotted points fall inside the displayed clock range"
  ),
  passed = c(
    !anyDuplicated(sleep_diary_valid_day_keys[c("site", "Id", "waking_date")]),
    setequal(
      sleep_diary_nights$participant,
      paste(sleep_diary_primary$site, sleep_diary_primary$Id, sep = "::")
    ),
    all(
      abs(
        as.numeric(sleep_diary_nights$sleep_onset_wall) -
          as.numeric(sleep_diary_nights$sleepprep_wall) -
          sleep_diary_nights$sleepdelay * 60
      ) <
        1e-6
    ),
    all(
      sleep_diary_nights$sleep_period_hours > 0 &
        sleep_diary_nights$sleep_period_hours < 24
    ),
    as.numeric(as.POSIXct("2026-11-07 23:00:00", tz = "UTC") + 8 * 3600 / 2) ==
      as.numeric(as.POSIXct("2026-11-08 03:00:00", tz = "UTC")),
    isTRUE(all.equal((((c(23.5, 0.5) + 12) %% 24) - 12), c(-0.5, 0.5))),
    all(
      sleep_diary_valid$midsleep_unwrapped_hours > -12 &
        sleep_diary_valid$midsleep_unwrapped_hours < 12
    ) &&
      diff(range(sleep_diary_valid$midsleep_unwrapped_hours)) < 12,
    identical(
      as.character(sleep_diary_nights$daytype),
      ifelse(
        as.character(sleep_diary_nights$daytype2) == "a work day",
        "Work",
        "Free"
      )
    ),
    all(
      sleep_diary_participant_sd$nights[is.finite(
        sleep_diary_participant_sd$sd_minutes
      )] >
        2
    ),
    all(is.na(sleep_diary_participant_sd$sd_minutes[
      sleep_diary_participant_sd$nights <= 2
    ])),
    all(
      abs(sleep_diary_dst_audit$elapsed_minus_wall_minutes[
        !sleep_diary_dst_audit$dst_crossing
      ]) <
        1e-6
    ),
    all(is.finite(as.matrix(select(
      sleep_diary_correlation_table,
      ends_with("lower"),
      ends_with("upper")
    )))),
    all(
      sleep_diary_sd_summary$minimum_nights >= 3 &
        sleep_diary_sd_summary$maximum_nights <=
          max(sleep_diary_participant_sd$nights)
    ),
    !anyDuplicated(sleep_diary_scatter_means_data[c("participant", "daytype")]),
    sum(sleep_diary_scatter_means_data$nights) == nrow(sleep_diary_nights),
    all(vapply(
      c("Work", "Free"),
      function(context) {
        means <- filter(sleep_diary_scatter_means_data, daytype == context)
        expected <- filter(
          sleep_diary_correlation_table,
          scope == "Primary exposure participants",
          .data$context == .env$context
        )
        nrow(means) == expected$participants &&
          isTRUE(all.equal(
            stats::cor(
              means$msfsc_unwrapped_hours,
              means$midsleep_unwrapped_hours
            ),
            expected$participant_mean_r
          ))
      },
      logical(1)
    )),
    all(
      is.finite(sleep_diary_msfsc_participants$diary_msfsc_hours) ==
        sleep_diary_msfsc_participants$eligible
    ),
    all(vapply(
      seq_len(nrow(sleep_diary_msfsc_eligible)),
      function(i) {
        person <- sleep_diary_msfsc_eligible$participant[i]
        expected <- mean(sleep_diary_nights$sleep_period_hours[
          sleep_diary_nights$participant == person
        ])
        abs(
          expected - sleep_diary_msfsc_eligible$mean_sleep_duration_hours[i]
        ) <
          1e-10
      },
      logical(1)
    )),
    all(
      sleep_diary_msfsc_eligible$correction_hours[
        sleep_diary_msfsc_eligible$sd_free_hours <=
          sleep_diary_msfsc_eligible$sd_work_hours
      ] <
        1e-10
    ),
    all(
      sleep_diary_msfsc_eligible$correction_hours >= 0 &
        sleep_diary_msfsc_eligible$diary_msfsc_hours <=
          sleep_diary_msfsc_eligible$msf_hours
    ),
    all(is.finite(as.matrix(select(
      sleep_diary_msfsc_summary,
      r,
      lower,
      upper
    )))),
    !anyDuplicated(sleep_diary_msfsc_eligible$participant) &&
      sum(sleep_diary_scatter_data$panel == "C  Diary-derived MSFsc") ==
        nrow(sleep_diary_msfsc_eligible),
    sum(sleep_diary_scatter_histogram_data$participants) ==
      nrow(sleep_diary_scatter_means_data),
    all(
      sleep_diary_scatter_data$msfsc_unwrapped_hours >= 0 &
        sleep_diary_scatter_data$msfsc_unwrapped_hours <= 12 &
        sleep_diary_scatter_data$midsleep_unwrapped_hours >= 0 &
        sleep_diary_scatter_data$midsleep_unwrapped_hours <= 12
    )
  )
)
stopifnot(all(sleep_diary_checks$passed))
sleep_diary_checks |> count(passed, name = "checks")

## ---- sleep-diary-export
for (directory in c(
  "results/tables",
  "results/data",
  "results/figures",
  "results/provenance",
  "data/intermediate"
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
sleep_diary_outputs <- list(
  sample_flow = sleep_diary_flow,
  correlations = sleep_diary_correlation_table,
  participant_sd = sleep_diary_participant_sd,
  sd_summary = sleep_diary_sd_summary,
  msfsc_participants = sleep_diary_msfsc_participants,
  msfsc_flow = sleep_diary_msfsc_flow,
  msfsc_summary = sleep_diary_msfsc_summary,
  exclusions = sleep_diary_records |>
    count(in_primary_sample, exclusion, name = "diary_entries"),
  dst_sensitivity = sleep_diary_dst_sensitivity
)
for (name in names(sleep_diary_outputs)) {
  write.csv(
    sleep_diary_outputs[[name]],
    paste0("results/tables/sleep_diary_", name, ".csv"),
    row.names = FALSE,
    na = ""
  )
}
write.csv(
  sleep_diary_records |>
    select(
      site,
      Id,
      source_row,
      participant,
      waking_date,
      daytype,
      msfsc_hours,
      midsleep_hours,
      midsleep_unwrapped_hours,
      difference_minutes,
      sleep_period_hours,
      in_primary_sample,
      waking_date_in_light_analysis,
      dst_crossing,
      exclusion
    ),
  "results/data/sleep_diary_nights.csv",
  row.names = FALSE,
  na = ""
)
write.csv(
  sleep_diary_scatter_data,
  "data/intermediate/sleep_diary_scatter_data.csv",
  row.names = FALSE,
  na = ""
)
write.csv(
  sleep_diary_scatter_correlations_data,
  "data/intermediate/sleep_diary_scatter_correlations_data.csv",
  row.names = FALSE
)
write.csv(
  sleep_diary_scatter_histogram_data,
  "data/intermediate/sleep_diary_scatter_histogram_data.csv",
  row.names = FALSE
)
ggsave(
  "results/figures/sleep_diary_scatter.png",
  sleep_diary_scatter,
  width = 13,
  height = 5,
  dpi = 220,
  bg = "white"
)
ggsave(
  "results/figures/sleep_diary_scatter.pdf",
  sleep_diary_scatter,
  width = 13,
  height = 5,
  bg = "white"
)
write.csv(
  sleep_diary_checks,
  "results/provenance/sleep_diary_checks.csv",
  row.names = FALSE
)
write.csv(
  sleep_diary_dst_audit |>
    filter(dst_crossing) |>
    select(
      site,
      Id,
      waking_date,
      daytype,
      midsleep_hours,
      elapsed_midpoint_hours,
      elapsed_minus_wall_minutes
    ),
  "results/provenance/sleep_diary_dst_audit.csv",
  row.names = FALSE
)
sleep_diary_input_files <- c(
  sleep_diary_manifest$path,
  "results/data/prepared.rds",
  "R/sleep_diary_alignment.R"
)
write.csv(
  tibble::tibble(
    path = sleep_diary_input_files,
    md5 = unname(tools::md5sum(sleep_diary_input_files))
  ),
  "results/provenance/sleep_diary_inputs.csv",
  row.names = FALSE
)
writeLines(
  c(
    "Command: Rscript --vanilla R/sleep_diary_alignment.R",
    "Also executed in analysis.qmd via knitr::read_chunk at the corresponding analytical stages.",
    "Seed: 20261106; 2000 participant-cluster bootstrap replicates per scope/context.",
    "Diary MSFsc-equivalent seed: 20261107; 2000 participant bootstrap replicates on the common eligible sample.",
    "Diary correction uses mean sleep-period duration weighted by each participant's observed work/free night counts; correction is zero if free-day duration is not longer.",
    "No diary alarm-use field is available. This is an MSFsc equivalent, not a full MCTQ rescore.",
    "Upstream daytype2 field label: Today is...; wake label: final awakening today.",
    "Sleep onset field verified as sleepprep_wall plus sleepdelay minutes.",
    "Wall-clock timestamps encode local labels in UTC; UTC fields retain elapsed-time coordinates.",
    "Lifecycle references checked 2026-11-06: dplyr.tidyverse.org/reference/summarise.html; dplyr.tidyverse.org/reference/mutate-joins.html; ggplot2.tidyverse.org/reference/geom_point.html.",
    "Formula checked 2026-11-06: docs.ropensci.org/mctq/reference/msf_sc.html; docs.ropensci.org/mctq/reference/sd_week.html. Observed night counts replace the usual weekly weights.",
    capture.output(sessionInfo())
  ),
  "results/provenance/sleep_diary_runtime.txt"
)
stopifnot(
  nrow(read.csv("data/intermediate/sleep_diary_scatter_data.csv")) ==
    nrow(sleep_diary_scatter_data)
)
