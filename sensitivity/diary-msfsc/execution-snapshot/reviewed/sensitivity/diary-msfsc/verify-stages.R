# Run from one isolated scenario directory after the primary replay finishes.
# Reuse the established numerical acceptance tests, adapting only their fixed
# original-sample assertion to the scenario's matched baseline.
source("R/setup.R")
for (module in c(
  "data",
  "models",
  "inference",
  "diagnostics",
  "validation",
  "noncyclic_sz"
))
  source(paste0("R/", module, ".R"))
prepared <- readRDS("results/data/prepared.rds")
checks <- list()
add <- function(name, passed, detail = "") {
  checks[[length(checks) + 1L]] <<- tibble(
    check = name,
    passed = isTRUE(passed),
    detail
  )
}
for (axis in c("clock", "phase")) {
  baseline <- readRDS(paste0("results/models/primary_", axis, ".rds"))
  noncyclic <- readRDS(paste0("results/models/noncyclic_sz_", axis, ".rds"))
  original_save_table <- save_table
  save_table <- function(...) invisible(NULL)
  basis_checks <- noncyclic_sz_acceptance(noncyclic, baseline)
  save_table <- original_save_table
  sample_row <- which(grepl("^Sample: 140 participants", basis_checks$check))
  stopifnot(length(sample_row) == 1L)
  basis_checks$check[
    sample_row
  ] <- "Identical sample counts to matched cyclic baseline"
  basis_checks$error[sample_row] <- sum(abs(
    unlist(count_sample(noncyclic$data)) - unlist(count_sample(baseline$data))
  ))
  basis_checks$passed[sample_row] <- basis_checks$error[sample_row] == 0
  save_table(
    basis_checks,
    paste0("noncyclic_sz_", axis, "_scenario_acceptance"),
    "results/provenance"
  )
  add(
    paste(
      axis,
      "noncyclic basis definitions, constraints and retained cyclic terms"
    ),
    all(basis_checks$passed),
    paste(nrow(basis_checks), "numerical checks")
  )
  shapley <- read.csv(paste0("results/tables/shapley_", axis, ".csv"))
  subsets <- read.csv(paste0("results/tables/shapley_subsets_", axis, ".csv"))
  add(
    paste(axis, "all 128 Shapley subsets retain the full-model AR coefficient"),
    identical(subsets$mask, 0:127) &&
      all(subsets$rho == subsets$rho[1]) &&
      max(abs(subsets$rho - baseline$rho)) < 1e-12,
    "CSV/RDS numerical comparison tolerance: 1e-12"
  )
  add(
    paste(axis, "Shapley allocations sum to the full-minus-baseline R-squared"),
    abs(
      sum(shapley$allocation) - (shapley$full_R2[1] - shapley$baseline_R2[1])
    ) <
      1e-9
  )
  rm(baseline, noncyclic)
  gc()
}
for (label in c("primary", "chest"))
  for (axis in c("clock", "phase")) {
    inf <- readRDS(paste0("results/data/", label, "_", axis, "_inference.rds"))
    add(
      paste(label, axis, "registered simultaneous coverage and support"),
      identical(inf$axis_coverage, cfg$curve_coverage) &&
        abs(
          inf$critical -
            unname(quantile(inf$simulated_maxima, cfg$curve_coverage, type = 8))
        ) <
          1e-12 &&
        nrow(inf$registry$W) == 20L &&
        all(abs(rowSums(inf$registry$W)) < 1e-10) &&
        !any(inf$contrasts$significant & !inf$contrasts$supported)
    )
  }
for (stage in c(
  "primary",
  "validation",
  "attribution",
  "sensitivity",
  "outputs",
  "meq",
  "daylight",
  "daylight-validation",
  "noncyclic",
  "prc"
)) {
  path <- paste0("results/provenance/", stage, "_run.csv")
  add(
    paste(stage, "completed with matching input signature"),
    file.exists(path) &&
      identical(read.csv(path)$prepared_signature, prepared$signature)
  )
}
checks <- bind_rows(checks)
save_table(checks, "scenario_stage_checks", "results/provenance")
capture.output(
  sessionInfo(),
  file = "results/provenance/scenario_stage_checks_session.txt"
)
stopifnot(all(checks$passed))
log_step("All", nrow(checks), "scenario-stage checks passed")
