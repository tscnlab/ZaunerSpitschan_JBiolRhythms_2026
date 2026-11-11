# Run from the chronotype_analysis project root:
# Rscript --vanilla sensitivity/diary-msfsc/prepare.R
source("R/setup.R")
source("R/data.R")
original_root <- normalizePath(getwd())
sensitivity_root <- file.path(original_root, "sensitivity/diary-msfsc")
dir.create(sensitivity_root, recursive = TRUE, showWarnings = FALSE)

## ---- dmsfsc-input-provenance
input_files <- unique(c(
  list.files("R", pattern = "[.]R$", full.names = TRUE),
  ".Rprofile",
  "renv.lock",
  "input-data/manifest.csv",
  list.files("input-data", recursive = TRUE, full.names = TRUE),
  "results/data/prepared.rds",
  Sys.glob("results/models/primary_*.rds"),
  Sys.glob("results/models/chest_*.rds"),
  Sys.glob("results/data/primary_*_inference.rds"),
  Sys.glob("results/tables/primary_*_significant_windows.csv")
))
input_files <- input_files[file.exists(input_files) & !dir.exists(input_files)]
hashes <- tibble(
  path = normalizePath(input_files),
  md5 = unname(tools::md5sum(input_files))
)
hash_file <- file.path(sensitivity_root, "original-input-hashes.csv")
if (file.exists(hash_file)) {
  previous <- read.csv(hash_file)
  stopifnot(
    identical(previous$path, hashes$path),
    identical(previous$md5, hashes$md5)
  )
} else write.csv(hashes, hash_file, row.names = FALSE)
snapshot <- file.path(sensitivity_root, "source-snapshot")
dir.create(snapshot, showWarnings = FALSE)
for (file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  target <- file.path(snapshot, basename(file))
  if (!file.exists(target)) stopifnot(file.copy(file, target))
  stopifnot(unname(tools::md5sum(file)) == unname(tools::md5sum(target)))
}

## ---- dmsfsc-derivation
paths <- input_paths()
chronotype <- readRDS(paths[["chronotype"]]) |>
  transmute(site, Id, questionnaire_msfsc_hours = as.numeric(msf_sc) / 3600)
diary <- readRDS(paths[["sleepdiaries"]]) |>
  left_join(chronotype, by = c("site", "Id"), relationship = "many-to-one") |>
  mutate(
    participant_key = paste(site, Id, sep = "::"),
    waking_date = as.Date(wake_wall, tz = "UTC"),
    daytype = factor(
      as.character(daytype2),
      levels = c("a work day", "a free day"),
      labels = c("Work", "Free")
    ),
    duration_hours = as.numeric(difftime(
      wake_wall,
      sleep_onset_wall,
      units = "hours"
    )),
    midpoint_wall = sleep_onset_wall + duration_hours * 3600 / 2,
    midsleep_hours = (as.numeric(midpoint_wall) / 3600) %% 24,
    midsleep_unwrapped = ((midsleep_hours + 12) %% 24) - 12,
    eligible_night = is.finite(questionnaire_msfsc_hours) &
      !is.na(sleep_onset_wall) &
      !is.na(wake_wall) &
      sleep_interval_analysis_eligible &
      !sleep_interval_quarantined &
      duration_hours > 0 &
      duration_hours < 24 &
      !is.na(daytype)
  )
valid <- diary |> filter(eligible_night)
stopifnot(!anyDuplicated(valid[c("participant_key", "waking_date")]))
derived <- valid |>
  summarise(
    nights = n(),
    midsleep = mean(midsleep_unwrapped),
    sleep_duration = mean(duration_hours),
    .by = c(site, Id, participant_key, daytype)
  ) |>
  pivot_wider(
    names_from = daytype,
    values_from = c(nights, midsleep, sleep_duration)
  ) |>
  mutate(
    eligible = !is.na(nights_Work) &
      !is.na(nights_Free) &
      nights_Work >= 1L &
      nights_Free >= 1L,
    weighted_sleep_duration = (nights_Work *
      sleep_duration_Work +
      nights_Free * sleep_duration_Free) /
      (nights_Work + nights_Free),
    correction_hours = pmax(0, sleep_duration_Free - weighted_sleep_duration) /
      2,
    diary_msfsc_hours = midsleep_Free - correction_hours
  ) |>
  left_join(chronotype, by = c("site", "Id"), relationship = "one-to-one") |>
  mutate(
    offset_minutes = wrap12(diary_msfsc_hours - questionnaire_msfsc_hours) * 60
  )
write.csv(
  diary,
  file.path(sensitivity_root, "diary-night-registry.csv"),
  row.names = FALSE
)
write.csv(
  derived,
  file.path(sensitivity_root, "derived-msfsc-participants.csv"),
  row.names = FALSE
)
parent_derivation <- read.csv(
  "results/tables/sleep_diary_msfsc_participants.csv"
) |>
  filter(eligible) |>
  select(site, Id, diary_msfsc_hours)
replication <- inner_join(
  filter(derived, eligible),
  parent_derivation,
  by = c("site", "Id"),
  suffix = c("", "_reported"),
  relationship = "one-to-one"
)
stopifnot(
  nrow(replication) == nrow(parent_derivation),
  max(abs(
    replication$diary_msfsc_hours - replication$diary_msfsc_hours_reported
  )) <
    1e-10
)

## ---- dmsfsc-scenario-data
original <- readRDS("results/data/prepared.rds")
eligible <- filter(derived, eligible)
scenarios <- c("matched-msfsc", "dmsfsc")
for (scenario in scenarios) {
  scenario_root <- file.path(sensitivity_root, scenario)
  dir.create(
    file.path(scenario_root, "R"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  for (folder in c(
    "data",
    "tables",
    "models",
    "figures",
    "diagnostics",
    "subsets",
    "validation",
    "provenance",
    "logs"
  ))
    dir.create(
      file.path(scenario_root, "results", folder),
      recursive = TRUE,
      showWarnings = FALSE
    )
  for (file in list.files(snapshot, full.names = TRUE))
    if (!file.exists(file.path(scenario_root, "R", basename(file))))
      stopifnot(file.copy(file, file.path(scenario_root, "R", basename(file))))
  if (!file.exists(file.path(scenario_root, "input-data")))
    stopifnot(file.symlink(
      file.path(original_root, "input-data"),
      file.path(scenario_root, "input-data")
    ))
  stopifnot(file.copy(
    ".Rprofile",
    file.path(scenario_root, ".Rprofile"),
    overwrite = TRUE
  ))
  frames <- lapply(original$frames, function(x) {
    z <- x |>
      inner_join(
        select(eligible, participant_key, diary_msfsc_hours),
        by = "participant_key",
        relationship = "many-to-one"
      ) |>
      mutate(
        questionnaire_msfsc_hours = msf_hours,
        msf_hours = if (scenario == "dmsfsc") diary_msfsc_hours else msf_hours,
        phase = wrap24(clock - msf_hours),
        phase_display = wrap12(clock - msf_hours)
      ) |>
      assign_sequences()
    stopifnot(
      all(is.finite(z$msf_hours)),
      !anyDuplicated(z[c("participant_day_key", "clock_bin")])
    )
    z
  })
  cutoffs <- quartile_cutoffs(frames$glasses)
  frames <- lapply(frames, apply_quartiles, cutoffs = cutoffs)
  registry <- bind_rows(frames, .id = "placement") |>
    select(
      all_of(names(original$registry)),
      questionnaire_msfsc_hours,
      diary_msfsc_hours
    )
  flow <- bind_rows(lapply(names(frames), function(placement) {
    bind_rows(
      mutate(
        count_sample(original$frames[[placement]]),
        stage = "Original analysis sample"
      ),
      mutate(
        count_sample(frames[[placement]]),
        stage = "Both diary day types available"
      )
    ) |>
      mutate(placement = placement, .before = 1)
  }))
  participants <- distinct(
    frames$glasses,
    participant_key,
    quartile,
    direction,
    extremity,
    msf_hours,
    meq,
    meq_type
  )
  prepared <- list(
    frames = frames,
    registry = registry,
    cutoffs = cutoffs,
    manifest = original$manifest,
    flow = flow,
    participants = participants,
    scenario = scenario,
    derivation = eligible,
    signature = digest::digest(list(scenario, frames, eligible))
  )
  saveRDS(
    prepared,
    file.path(scenario_root, "results/data/prepared.rds"),
    compress = FALSE
  )
  write.csv(
    registry,
    file.path(scenario_root, "results/data/sample_registry.csv"),
    row.names = FALSE
  )
  write.csv(
    flow,
    file.path(scenario_root, "results/tables/sample_flow.csv"),
    row.names = FALSE
  )
  write.csv(
    tibble(
      percentile = c(25, 50, 75),
      msf_hours = cutoffs,
      clock_time = clock_label(cutoffs)
    ),
    file.path(scenario_root, "results/tables/quartile_cutoffs.csv"),
    row.names = FALSE
  )
  write.csv(
    participants,
    file.path(scenario_root, "results/tables/participants.csv"),
    row.names = FALSE
  )
  write.csv(
    bind_rows(frames, .id = "placement") |>
      summarise(
        participants = n_distinct(participant_key),
        participant_days = n_distinct(participant_day_key),
        half_hours = n(),
        .by = c(placement, quartile, daytype)
      ),
    file.path(scenario_root, "results/tables/group_day_counts.csv"),
    row.names = FALSE
  )
}
matched <- readRDS(file.path(
  sensitivity_root,
  "matched-msfsc/results/data/prepared.rds"
))
changed <- readRDS(file.path(
  sensitivity_root,
  "dmsfsc/results/data/prepared.rds"
))
for (placement in names(matched$frames)) {
  unchanged <- c(
    "participant_key",
    "participant_day_key",
    "clock_bin",
    "clock",
    "response",
    "daytype",
    "photoperiod_hours",
    "AR_start",
    "source_utc_start",
    "true_elapsed_sequence_id"
  )
  stopifnot(identical(
    matched$frames[[placement]][unchanged],
    changed$frames[[placement]][unchanged]
  ))
}
overlap <- inner_join(
  select(
    matched$participants,
    participant_key,
    questionnaire_quartile = quartile
  ),
  select(changed$participants, participant_key, diary_quartile = quartile),
  by = "participant_key",
  relationship = "one-to-one"
)
write.csv(
  overlap,
  file.path(sensitivity_root, "quartile-reassignment-participants.csv"),
  row.names = FALSE
)
write.csv(
  count(overlap, questionnaire_quartile, diary_quartile),
  file.path(sensitivity_root, "quartile-reassignment.csv"),
  row.names = FALSE
)
write.csv(
  bind_rows(lapply(
    scenarios,
    function(scenario)
      mutate(
        read.csv(file.path(
          sensitivity_root,
          scenario,
          "results/tables/sample_flow.csv"
        )),
        scenario = scenario
      )
  )),
  file.path(sensitivity_root, "sample-flow.csv"),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(sensitivity_root, "preparation-session.txt")
)
cat(
  "Prepared matched and dMSFsc scenarios:",
  nrow(matched$participants),
  "near-eye participants.\n"
)
