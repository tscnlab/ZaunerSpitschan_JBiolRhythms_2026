# Refresh target-level counts after removing unused factor levels from export.
source("R/setup.R")
for (module in c("data", "models", "inference", "daylight"))
  source(paste0("R/", module, ".R"))
source("../adapters.R")
bundle <- readRDS("results/models/daylight_phase.rds")
for (specification in c("fixed", "quartiles", "p10_p90")) {
  prefix <- if (specification == "fixed") "daylight" else
    paste0("daylight_", specification)
  old <- read.csv(paste0("results/tables/", prefix, "_scenario_samples.csv")) |>
    filter(quartile %in% c("Q1", "Q4"))
  inference_path <- paste0("results/data/", prefix, "_inference.rds")
  before <- tools::md5sum(inference_path)
  design <- diary_daylight_design(bundle$data, specification, prefix)
  current <- read.csv(paste0(
    "results/tables/",
    prefix,
    "_scenario_samples.csv"
  ))
  old <- arrange(old, scenario, quartile, daytype)
  current <- arrange(current, scenario, quartile, daytype)
  rownames(old) <- rownames(current) <- NULL
  stopifnot(
    nrow(current) == 12L,
    identical(old, current),
    identical(before, tools::md5sum(inference_path))
  )
}
cat("Target support counts verified; inference files unchanged.\n")
