# Fit phase and daylength models with diary timing and fixed questionnaire groups.
# Run from dmsfsc-fixed-groups: Rscript --vanilla -e 'source("../fixed-groups.R")' primary
source("R/setup.R")
for (module in c("data", "models", "inference", "diagnostics", "daylight")) {
  source(paste0("R/", module, ".R"))
}
source("../adapters.R")
stage <- match.arg(
  commandArgs(trailingOnly = TRUE)[1],
  c("primary", "daylight")
)
prepared <- readRDS("results/data/prepared.rds")
started <- Sys.time()

## ---- fixed-groups-classification
# The shared preprocessing normally reassigns groups from msf_hours. Use the
# questionnaire value only for classification, then restore diary timing.
original_apply_quartiles <- apply_quartiles
apply_quartiles <- function(x, cutoffs) {
  stopifnot(all(is.finite(x$questionnaire_msfsc_hours)))
  diary_timing <- x$msf_hours
  x$msf_hours <- x$questionnaire_msfsc_hours
  x <- original_apply_quartiles(x, cutoffs)
  x$msf_hours <- diary_timing
  x
}
raw <- prepared$frames$glasses
prep <- fit_preprocess(raw, prepared$cutoffs)
model_data <- apply_preprocess(raw, prep, "phase")
reference <- readRDS("../matched-msfsc/results/models/primary_phase.rds")$data
diary_data <- readRDS("../dmsfsc/results/models/primary_phase.rds")$data
fixed_columns <- setdiff(
  names(reference),
  c("msf_hours", "phase", "phase_display", "t")
)
stopifnot(
  identical(model_data[fixed_columns], reference[fixed_columns]),
  identical(
    model_data[c("t", "phase", "phase_display", "msf_hours")],
    diary_data[c("t", "phase", "phase_display", "msf_hours")]
  )
)

## ---- fixed-groups-primary
if (stage == "primary") {
  bundle <- fit_pair(raw, prep, "phase", "primary")
  stopifnot(identical(bundle$data$quartile, reference$quartile))
  inference <- infer_curves(bundle)
  inference$individual_draws <- NULL
  atomic_rds(inference, "results/data/primary_phase_inference.rds")
  save_table(
    interval_average_contrasts(inference, bundle),
    "primary_phase_interval_averages"
  )
  diagnose_model(bundle)
  save_table(model_summary(bundle), "primary_phase_model_fit")
}

## ---- fixed-groups-daylength
if (stage == "daylight") {
  bundle <- daylight_fit(raw, prep)
  stopifnot(identical(bundle$data$quartile, reference$quartile))
  basis_checks <- daylight_basis_checks(bundle)
  stopifnot(all(basis_checks$passed))
  save_table(basis_checks, "daylight_basis_checks", "results/provenance")
  diagnose_model(bundle)
  save_table(model_summary(bundle), "daylight_model_summary")
  windows <- read.csv(
    "../../../results/tables/primary_phase_significant_windows.csv"
  ) |>
    filter(contrast == "evening_morning_Extreme_Overall") |>
    mutate(window_source = "Original primary fixed")
  interval_results <- list()
  for (specification in c("fixed", "quartiles", "p10_p90")) {
    prefix <- if (specification == "fixed") "daylight" else
      paste0("daylight_", specification)
    design <- diary_daylight_design(bundle$data, specification, prefix)
    inference <- daylight_inference(bundle, design)
    daylight_windows(bundle, inference, design)
    interval_results[[specification]] <- diary_daylight_intervals(
      bundle,
      inference,
      design,
      windows
    ) |>
      mutate(target_specification = specification)
  }
  save_table(bind_rows(interval_results), "daylight_interval_comparison")
}

## ---- fixed-groups-run-record
inputs <- read.csv("results/provenance/input-hashes.csv")
stopifnot(identical(inputs$md5, unname(tools::md5sum(inputs$path))))
capture.output(
  sessionInfo(),
  file = paste0("results/provenance/", stage, "_session.txt")
)
save_table(
  tibble(
    stage,
    scenario = prepared$scenario,
    started = as.character(started),
    finished = as.character(Sys.time()),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    prepared_signature = prepared$signature,
    seed = cfg$seed,
    covariance_draws = cfg$draws,
    support_n = cfg$support_n,
    within_coordinate_coverage = cfg$curve_coverage
  ),
  paste0(stage, "_run"),
  "results/provenance"
)
log_step("Completed", prepared$scenario, stage)
