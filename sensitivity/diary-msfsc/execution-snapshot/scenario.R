# Invoke with a scenario directory as the working directory:
# Rscript --vanilla ../scenario.R primary
source("R/setup.R")
for (module in c(
  "data",
  "models",
  "inference",
  "diagnostics",
  "validation",
  "attribution",
  "plots",
  "sensitivity",
  "meq",
  "daylight",
  "phase_response",
  "emens_balance",
  "noncyclic_sz"
))
  source(paste0("R/", module, ".R"))
arguments <- commandArgs(trailingOnly = TRUE)
stage <- if (length(arguments)) arguments[1] else "primary"
prepared <- readRDS("results/data/prepared.rds")
scenario <- prepared$scenario
axis_title <- function(axis)
  if (axis == "clock") "Local clock time" else
    paste("Hours relative to", if (scenario == "dmsfsc") "dMSFsc" else "MSFsc")
started <- Sys.time()
source("../adapters.R")

## ---- dmsfsc-primary-fits
if (stage == "primary") {
  for (label in c("primary", "chest")) {
    raw <- prepared$frames[[if (label == "primary") "glasses" else "chest"]]
    prep <- fit_preprocess(raw, prepared$cutoffs)
    for (axis in c("clock", "phase")) {
      bundle <- fit_pair(raw, prep, axis, label)
      inference <- infer_curves(bundle)
      inference$individual_draws <- NULL
      atomic_rds(
        inference,
        paste0("results/data/", label, "_", axis, "_inference.rds")
      )
      intervals <- interval_average_contrasts(inference, bundle)
      save_table(intervals, paste0(label, "_", axis, "_interval_averages"))
      diagnose_model(bundle)
      save_table(model_summary(bundle), paste0(label, "_", axis, "_model_fit"))
      rm(bundle, inference)
      gc()
    }
  }
}

## ---- dmsfsc-validation
if (stage == "validation") run_validation(prepared)

## ---- dmsfsc-attribution
if (stage == "attribution") run_attribution(prepared)

## ---- dmsfsc-model-sensitivities
if (stage == "sensitivity") run_sensitivity(prepared)

## ---- dmsfsc-descriptive-outputs
if (stage == "outputs") {
  for (label in c("primary", "chest")) {
    curves <- contrasts <- list()
    for (axis in c("clock", "phase")) {
      bundle <- readRDS(paste0("results/models/", label, "_", axis, ".rds"))
      inference <- readRDS(paste0(
        "results/data/",
        label,
        "_",
        axis,
        "_inference.rds"
      ))
      curves[[axis]] <- mutate(inference$curves, axis = axis)
      contrasts[[axis]] <- mutate(inference$contrasts, axis = axis)
      photoperiod_figure(bundle)
      if (label == "primary") individual_curves(bundle, inference)
      rm(bundle, inference)
      gc()
    }
    for (context in c("Overall", "Work", "Free")) {
      save_plot(
        pattern_figure(bind_rows(curves), context, label),
        paste0(label, "_patterns_", tolower(context)),
        12,
        11
      )
      for (family in c("Evening minus morning", "Extreme minus moderate"))
        save_plot(
          difference_figure(bind_rows(contrasts), family, context, label),
          paste0(
            label,
            "_",
            gsub(" ", "_", tolower(family)),
            "_",
            tolower(context)
          ),
          12,
          11
        )
    }
    for (family in c("Free minus work", "Interactions"))
      save_plot(
        difference_figure(bind_rows(contrasts), family, placement = label),
        paste0(label, "_", gsub(" ", "_", tolower(family))),
        12,
        14
      )
  }
  matched_placement_check(prepared)
}

## ---- dmsfsc-meq
if (stage == "meq") {
  raw <- meq_prepare(prepared)
  prep <- fit_preprocess(raw, prepared$cutoffs)
  curves <- contrasts <- windows <- summaries <- checks <- list()
  for (classification in c("MEQ", "MSFsc"))
    for (axis in c("clock", "phase")) {
      key <- paste(classification, axis, sep = "_")
      bundle <- meq_fit(raw, prep, axis, classification)
      inference <- meq_inference(bundle)
      summaries[[key]] <- mutate(
        model_summary(bundle),
        classification = classification
      )
      curves[[key]] <- inference$curves
      contrasts[[key]] <- inference$contrasts
      windows[[key]] <- inference$windows
      diagnose_model(bundle)
      checks[[key]] <- tibble(
        classification,
        axis,
        observations_match = identical(
          bundle$data$participant_day_key,
          raw$participant_day_key
        ) &&
          identical(bundle$data$clock_bin, raw$clock_bin),
        support_valid = all(
          !inference$contrasts$candidate_interval |
            inference$contrasts$minimum_cell_participants >= cfg$support_n
        )
      )
      rm(bundle, inference)
      gc()
    }
  save_table(bind_rows(curves), "meq_all_curves")
  save_table(bind_rows(contrasts), "meq_all_contrasts")
  save_table(bind_rows(windows), "meq_candidate_windows")
  save_table(bind_rows(summaries), "meq_model_fit")
  checks <- bind_rows(checks)
  stopifnot(all(checks$observations_match), all(checks$support_valid))
  save_table(checks, "meq_acceptance", "results/provenance")
}

## ---- dmsfsc-daylength
if (stage == "daylight") {
  raw <- prepared$frames$glasses
  prep <- fit_preprocess(raw, prepared$cutoffs)
  bundle <- daylight_fit(raw, prep)
  checks <- daylight_basis_checks(bundle)
  stopifnot(all(checks$passed))
  diagnose_model(bundle)
  save_table(model_summary(bundle), "daylight_model_summary")
  own_windows <- read.csv(
    "results/tables/primary_phase_significant_windows.csv"
  ) |>
    filter(contrast == "evening_morning_Extreme_Overall")
  original_windows <- read.csv(
    "../../../results/tables/primary_phase_significant_windows.csv"
  ) |>
    filter(contrast == "evening_morning_Extreme_Overall")
  windows <- bind_rows(
    mutate(own_windows, window_source = "Scenario-selected"),
    mutate(original_windows, window_source = "Original primary fixed")
  )
  if (scenario == "dmsfsc") {
    matched_windows <- read.csv(
      "../matched-msfsc/results/tables/primary_phase_significant_windows.csv"
    ) |>
      filter(contrast == "evening_morning_Extreme_Overall")
    windows <- bind_rows(
      windows,
      mutate(matched_windows, window_source = "Matched questionnaire fixed")
    )
  }
  all_intervals <- list()
  for (specification in c("fixed", "quartiles", "p10_p90")) {
    prefix <- if (specification == "fixed") "daylight" else
      paste0("daylight_", specification)
    design <- diary_daylight_design(bundle$data, specification, prefix)
    inference <- daylight_inference(bundle, design)
    daylight_windows(bundle, inference, design)
    all_intervals[[specification]] <- diary_daylight_intervals(
      bundle,
      inference,
      design,
      windows
    ) |>
      mutate(target_specification = specification)
  }
  save_table(bind_rows(all_intervals), "daylight_interval_comparison")
  diary_solar_summary(bundle, windows)
}
if (stage == "daylight-validation") daylight_validation(prepared)

## ---- dmsfsc-noncyclic
if (stage == "noncyclic") {
  for (axis in c("clock", "phase")) {
    bundle <- noncyclic_sz_fit(axis)
    inference <- infer_curves(bundle)
    inference$windows <- noncyclic_sz_windows(inference, bundle)
    save_table(
      inference$windows,
      paste0("noncyclic_sz_", axis, "_significant_windows")
    )
    atomic_rds(
      inference,
      paste0("results/data/noncyclic_sz_", axis, "_inference.rds")
    )
    boundary <- noncyclic_sz_boundary(bundle)
    save_table(
      boundary$curves,
      paste0("noncyclic_sz_", axis, "_boundary_curves")
    )
    save_table(
      boundary$contrasts,
      paste0("noncyclic_sz_", axis, "_boundary_contrasts")
    )
    save_table(
      model_summary(bundle),
      paste0("noncyclic_sz_", axis, "_model_fit")
    )
    diagnose_model(bundle)
    rm(bundle, inference)
    gc()
  }
}

## ---- dmsfsc-phase-response
if (stage == "prc") {
  prc <- phase_response_read_prc()
  bundle <- readRDS("results/models/primary_phase.rds")
  inference <- readRDS("results/data/primary_phase_inference.rds")
  result <- phase_response_integrate(prc, bundle, inference)
  phase_response_export(result)
  emens_balance_export(emens_balance(bundle, inference))
}

capture.output(
  sessionInfo(),
  file = paste0("results/provenance/", stage, "_session.txt")
)
save_table(
  tibble(
    stage = stage,
    scenario = scenario,
    started = as.character(started),
    finished = as.character(Sys.time()),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    prepared_signature = prepared$signature,
    seed = cfg$seed,
    covariance_draws = cfg$draws,
    bootstrap_draws = cfg$bootstrap,
    support_n = cfg$support_n,
    within_coordinate_coverage = cfg$curve_coverage
  ),
  paste0(stage, "_run"),
  "results/provenance"
)
log_step("Completed", scenario, stage)
