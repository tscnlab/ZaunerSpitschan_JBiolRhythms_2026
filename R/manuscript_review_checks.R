# Descriptive checks requested during the November 2026 coauthor review.
# Existing model fits and inference are not changed.
source("R/setup.R")
source_demographics <- "../ZaunerEtAl_bioRxiv_2026/results/intermediate/model_data/normalized_inputs/demographics.rds"
snapshot <- "input-data/review-demographics.rds"
if (!file.exists(snapshot)) {
  if (!file.exists(source_demographics)) stop("Demographic source snapshot is unavailable")
  file.copy(source_demographics, snapshot)
}
demographics <- readRDS(snapshot)
participants <- readRDS("results/data/prepared.rds")$frames$glasses |>
  distinct(site, Id)
stopifnot(nrow(participants) == 140L,
  !anyDuplicated(demographics[c("site", "Id")]))
primary_demographics <- left_join(participants, demographics,
  by = c("site", "Id"), relationship = "one-to-one")
stopifnot(!anyNA(primary_demographics$age), !anyNA(primary_demographics$sex))
demographic_summary <- primary_demographics |>
  summarise(participants = n(), female = sum(sex == "Female"),
    male = sum(sex == "Male"), median_age = median(age),
    minimum_age = min(age), maximum_age = max(age))
write.csv(demographic_summary, "results/tables/primary_demographics.csv", row.names = FALSE)
reassignment <- read.csv("sensitivity/diary-msfsc/quartile-reassignment.csv")
transition_summary <- reassignment |>
  filter(questionnaire_quartile != diary_quartile) |>
  mutate(distance = abs(match(questionnaire_quartile, paste0("Q", 1:4)) -
    match(diary_quartile, paste0("Q", 1:4)))) |>
  arrange(desc(n), questionnaire_quartile, diary_quartile)
stopifnot(sum(reassignment$n) == 135L, sum(transition_summary$n) == 50L)
write.csv(transition_summary, "results/tables/diary_quartile_transitions.csv", row.names = FALSE)
inputs <- c(snapshot, "results/data/prepared.rds", "sensitivity/diary-msfsc/quartile-reassignment.csv")
write.csv(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))),
  "results/provenance/manuscript_review_inputs.csv", row.names = FALSE)
capture.output({
  cat("Command: Rscript --vanilla R/manuscript_review_checks.R\n")
  cat("Original demographics source:", source_demographics, "\n")
  print(demographic_summary)
  print(transition_summary)
  print(sessionInfo())
}, file = "results/provenance/manuscript_review_runtime.txt")
print(demographic_summary)
print(transition_summary)
