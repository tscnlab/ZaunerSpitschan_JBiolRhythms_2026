# Prepare diary-based phase alignment with questionnaire quartiles held fixed.
# Run from the chronotype_analysis project root.
source("R/setup.R")
root <- "sensitivity/diary-msfsc"
source(file.path(root, "source-snapshot/data.R"))
matched <- readRDS(file.path(root, "matched-msfsc/results/data/prepared.rds"))
diary <- readRDS(file.path(root, "dmsfsc/results/data/prepared.rds"))
scenario_root <- file.path(root, "dmsfsc-fixed-groups")

## ---- fixed-groups-preparation
group_columns <- c("quartile", "direction", "extremity", "group_label")
raw <- diary$frames$glasses
reference <- matched$frames$glasses
changing_columns <- c("msf_hours", "phase", "phase_display", group_columns)
unchanged_columns <- setdiff(names(reference), changing_columns)
stopifnot(
  identical(raw[unchanged_columns], reference[unchanged_columns]),
  !anyDuplicated(raw[c("participant_day_key", "clock_bin")])
)
raw[group_columns] <- reference[group_columns]
stopifnot(
  identical(raw[group_columns], reference[group_columns]),
  identical(raw$phase, diary$frames$glasses$phase),
  all(abs(raw$phase - wrap24(raw$clock - raw$diary_msfsc_hours)) < 1e-12),
  all(abs(raw$msf_hours - raw$diary_msfsc_hours) < 1e-12)
)
prepared <- list(
  frames = list(glasses = raw),
  cutoffs = matched$cutoffs,
  participants = distinct(
    raw,
    participant_key,
    quartile,
    direction,
    extremity,
    questionnaire_msfsc_hours,
    diary_msfsc_hours
  ),
  scenario = "dmsfsc-fixed-groups",
  grouping_metric = "Questionnaire MSFsc on the matched sample",
  phase_metric = "Diary-derived dMSFsc",
  signature = digest::digest(list(raw, matched$cutoffs, "fixed-groups-1.0"))
)

## ---- fixed-groups-preparation-exports
for (folder in c(
  "R",
  "results/data",
  "results/models",
  "results/tables",
  "results/diagnostics",
  "results/figures",
  "results/provenance",
  "results/logs"
)) {
  dir.create(
    file.path(scenario_root, folder),
    recursive = TRUE,
    showWarnings = FALSE
  )
}
producer_files <- list.files(
  file.path(root, "source-snapshot"),
  full.names = TRUE
)
stopifnot(all(file.copy(
  producer_files,
  file.path(scenario_root, "R"),
  overwrite = TRUE
)))
stopifnot(file.copy(".Rprofile", scenario_root, overwrite = TRUE))
saveRDS(
  prepared,
  file.path(scenario_root, "results/data/prepared.rds"),
  compress = FALSE
)
write.csv(
  prepared$participants,
  file.path(scenario_root, "results/tables/participants.csv"),
  row.names = FALSE
)
write.csv(
  count_sample(raw),
  file.path(scenario_root, "results/tables/sample.csv"),
  row.names = FALSE
)
input_files <- c(
  file.path(root, c("matched-msfsc", "dmsfsc"), "results/data/prepared.rds"),
  file.path(
    root,
    "source-snapshot",
    c(
      "setup.R",
      "data.R",
      "models.R",
      "inference.R",
      "daylight.R",
      "diagnostics.R"
    )
  ),
  file.path(root, c("adapters.R", "prepare-fixed-groups.R", "fixed-groups.R"))
)
write.csv(
  tibble(
    path = normalizePath(input_files),
    md5 = unname(tools::md5sum(input_files))
  ),
  file.path(scenario_root, "results/provenance/input-hashes.csv"),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(scenario_root, "results/provenance/preparation-session.txt")
)
print(count_sample(raw))
cat("Questionnaire membership retained for every observation.\n")
