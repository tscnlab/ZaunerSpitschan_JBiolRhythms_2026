# Compare the primary analysis with matched questionnaire and diary scenarios.
# Run from the original project root after all scenario stages have completed.
source("R/setup.R")
source("R/inference.R")
source("R/data.R")
source("R/models.R")
source("R/diagnostics.R")
root <- "sensitivity/diary-msfsc"
dir.create(file.path(root, "comparison"), showWarnings = FALSE)
dir.create(
  file.path(root, "data/intermediate"),
  recursive = TRUE,
  showWarnings = FALSE
)
paths <- c(
  "Original questionnaire (140)" = ".",
  "Matched questionnaire (135)" = file.path(root, "matched-msfsc"),
  "Diary dMSFsc (135)" = file.path(root, "dmsfsc")
)
# Display notation only; stored scenario identifiers remain unchanged.
scenario_labels <- function(x) {
  lapply(x, function(label) {
    if (grepl("dMSFsc", label, fixed = TRUE))
      bquote("Diary " * MSF["sc(diary)"] * " (135)") else label
  })
}
comparison_table <- function(data, name) {
  write.csv(
    data,
    file.path(root, "comparison", paste0(name, ".csv")),
    row.names = FALSE
  )
  invisible(data)
}
read_scenarios <- function(file) {
  bind_rows(lapply(names(paths), function(scenario) {
    path <- file.path(paths[[scenario]], "results", file)
    if (!file.exists(path)) stop("Required scenario result missing: ", path)
    mutate(read.csv(path), scenario = scenario)
  }))
}

## ---- dmsfsc-comparison-sample
derived <- read.csv(file.path(root, "derived-msfsc-participants.csv"))
matched <- readRDS(file.path(paths[[2]], "results/data/prepared.rds"))
changed <- readRDS(file.path(paths[[3]], "results/data/prepared.rds"))
offsets <- derived |>
  filter(eligible, participant_key %in% matched$participants$participant_key)
comparison_table(
  tibble(
    participants = nrow(offsets),
    correlation = cor(
      offsets$questionnaire_msfsc_hours,
      offsets$diary_msfsc_hours
    ),
    median_offset_minutes = median(offsets$offset_minutes),
    lower_quartile_offset_minutes = quantile(offsets$offset_minutes, .25),
    upper_quartile_offset_minutes = quantile(offsets$offset_minutes, .75),
    median_absolute_offset_minutes = median(abs(offsets$offset_minutes))
  ),
  "metric_comparison"
)
comparison_table(
  bind_rows(lapply(names(paths), function(scenario) {
    prep <- readRDS(file.path(paths[[scenario]], "results/data/prepared.rds"))
    mutate(
      count_sample(prep$frames$glasses),
      scenario = scenario,
      cutoff25 = prep$cutoffs[1],
      cutoff50 = prep$cutoffs[2],
      cutoff75 = prep$cutoffs[3]
    )
  })),
  "samples_and_cutoffs"
)

## ---- dmsfsc-comparison-primary
curves <- contrasts <- windows <- fixed_intervals <- fits <- list()
original_windows <- read.csv(
  "results/tables/primary_phase_significant_windows.csv"
) |>
  filter(contrast == "evening_morning_Extreme_Overall") |>
  mutate(window_source = "Original primary fixed")
matched_windows <- read.csv(file.path(
  paths[[2]],
  "results/tables/primary_phase_significant_windows.csv"
)) |>
  filter(contrast == "evening_morning_Extreme_Overall") |>
  mutate(window_source = "Matched questionnaire fixed")
for (scenario in names(paths))
  for (axis in c("clock", "phase")) {
    bundle <- readRDS(file.path(
      paths[[scenario]],
      paste0("results/models/primary_", axis, ".rds")
    ))
    inference <- readRDS(file.path(
      paths[[scenario]],
      paste0("results/data/primary_", axis, "_inference.rds")
    ))
    key <- paste(scenario, axis)
    curves[[key]] <- mutate(inference$curves, scenario = scenario, axis = axis)
    contrasts[[key]] <- mutate(
      inference$contrasts,
      scenario = scenario,
      axis = axis
    )
    windows[[key]] <- mutate(
      inference$windows,
      scenario = scenario,
      axis = axis
    )
    fits[[key]] <- mutate(model_summary(bundle), scenario = scenario)
    if (axis == "phase") {
      selected <- bind_rows(
        original_windows,
        matched_windows,
        filter(
          inference$windows,
          contrast == "evening_morning_Extreme_Overall"
        ) |>
          mutate(window_source = "Scenario-selected")
      )
      cv <- covariance_info(bundle$fit)
      for (i in seq_len(nrow(selected)))
        for (context in c("Overall", "Work", "Free")) {
          id <- paste0("evening_morning_Extreme_", context)
          weights <- inference$registry$W[id, ]
          D <- Reduce(
            `+`,
            Map(function(X, w) X * w, inference$cell_design$X, weights)
          )
          ix <- which(
            wrap24(inference$cell_design$x - selected$start[i]) <
              selected$duration_hours[i] - 1e-10
          )
          estimate <- linear_summary(
            matrix(colMeans(D[ix, , drop = FALSE]), nrow = 1),
            coef(bundle$fit),
            cv$V
          )
          local_support <- filter(
            inference$contrasts,
            contrast == id,
            x %in% inference$cell_design$x[ix]
          )
          fixed_intervals[[length(fixed_intervals) + 1L]] <- bind_cols(
            selected[i, c("start", "end", "duration_hours", "window_source")],
            estimate
          ) |>
            mutate(
              scenario = scenario,
              context = context,
              phase_label = ifelse(end < 0, "Pre-MSFsc", "Post-MSFsc"),
              crosses_boundary = start >= end,
              ratio = 10^estimate,
              ratio_lower = 10^lower,
              ratio_upper = 10^upper,
              supported = all(local_support$supported),
              minimum_cell_participants = min(
                local_support$minimum_cell_participants
              ),
              interval_method = "Mean log contrast; selected window fixed; unadjusted 95% Gaussian interval"
            )
        }
    }
    rm(bundle, inference)
    gc()
  }
curves <- bind_rows(curves)
contrasts <- bind_rows(contrasts)
windows <- bind_rows(windows)
comparison_table(curves, "primary_curves")
comparison_table(contrasts, "primary_contrasts")
comparison_table(windows, "primary_significant_windows")
comparison_table(bind_rows(fixed_intervals), "phase_interval_averages")
comparison_table(bind_rows(fits), "model_fit")
comparison_table(
  bind_rows(lapply(
    c("clock", "phase"),
    function(axis)
      read_scenarios(paste0("diagnostics/primary_", axis, "_verdict.csv"))
  )),
  "primary_diagnostics"
)

## ---- dmsfsc-comparison-daylength
daylength <- bind_rows(lapply(names(paths)[-1], function(scenario) {
  mutate(
    read.csv(file.path(
      paths[[scenario]],
      "results/tables/daylight_interval_comparison.csv"
    )),
    scenario = scenario
  )
}))
comparison_table(daylength, "daylength_intervals")
comparison_table(
  read_scenarios("tables/daylight_validation_scores.csv"),
  "daylength_validation_scores"
)
comparison_table(
  read_scenarios("tables/daylight_validation_paired_differences.csv"),
  "daylength_validation_differences"
)
comparison_table(
  read_scenarios("tables/validation_scores.csv"),
  "validation_scores"
)

## ---- dmsfsc-paired-validation-comparison
matched_predictions <- read.csv(file.path(
  paths[[2]],
  "results/validation/held_out_predictions.csv"
))
diary_predictions <- read.csv(file.path(
  paths[[3]],
  "results/validation/held_out_predictions.csv"
))
pair_keys <- c(
  "target",
  "axis",
  "split",
  "site",
  "participant_key",
  "participant_day_key",
  "local_date",
  "clock_bin"
)
matched_predictions <- arrange(matched_predictions, across(all_of(pair_keys)))
diary_predictions <- arrange(diary_predictions, across(all_of(pair_keys)))
stopifnot(
  identical(matched_predictions[pair_keys], diary_predictions[pair_keys]),
  isTRUE(all.equal(matched_predictions$observed, diary_predictions$observed)),
  isTRUE(all.equal(matched_predictions$baseline, diary_predictions$baseline))
)
paired_predictions <- bind_rows(
  mutate(matched_predictions, measure = "Questionnaire"),
  mutate(diary_predictions, measure = "Diary")
)
person_scores <- paired_predictions |>
  summarise(
    mse = mean((predicted - observed)^2),
    null_mse = mean((baseline - observed)^2),
    .by = c(target, axis, measure, site, participant_key)
  )
paired_results <- list()
set.seed(cfg$seed + 102L)
for (axis in c("clock", "phase"))
  for (target in unique(person_scores$target)) {
    a <- filter(
      person_scores,
      .data$axis == .env$axis,
      .data$target == .env$target,
      measure == "Questionnaire"
    ) |>
      arrange(site, participant_key)
    b <- filter(
      person_scores,
      .data$axis == .env$axis,
      .data$target == .env$target,
      measure == "Diary"
    ) |>
      arrange(site, participant_key)
    stopifnot(
      identical(a$participant_key, b$participant_key),
      isTRUE(all.equal(a$null_mse, b$null_mse))
    )
    score <- function(x, ix)
      c(
        RMSE = sqrt(mean(x$mse[ix])),
        predictive_R2 = 1 - mean(x$mse[ix]) / mean(x$null_mse[ix])
      )
    indices <- seq_len(nrow(a))
    difference <- score(b, indices) - score(a, indices)
    strata <- split(indices, a$site)
    bootstrap <- replicate(cfg$bootstrap, {
      ix <- unlist(lapply(
        strata,
        function(i) i[sample.int(length(i), length(i), replace = TRUE)]
      ))
      score(b, ix) - score(a, ix)
    })
    ci <- apply(bootstrap, 1, quantile, c(.025, .975), type = 8)
    paired_results[[paste(axis, target)]] <- tibble(
      axis,
      target,
      metric = names(difference),
      diary_minus_questionnaire = as.numeric(difference),
      lower = ci[1, ],
      upper = ci[2, ],
      participants = nrow(a),
      bootstrap_replicates = cfg$bootstrap
    )
  }
comparison_table(bind_rows(paired_results), "paired_metric_validation")

## ---- dmsfsc-comparison-prc
prc <- read_scenarios("tables/phase_response_contrast_report.csv")
emens <- read_scenarios("tables/emens_balance_contrasts.csv")
comparison_table(prc, "prc_contrasts")
comparison_table(
  read_scenarios("tables/phase_response_contrasts.csv"),
  "prc_alignment_contrasts"
)
comparison_table(emens, "emens_contrasts")
comparison_table(
  read_scenarios("tables/phase_response_overall_report.csv"),
  "prc_overall_indices"
)
comparison_table(
  read_scenarios("tables/emens_balance_summary.csv"),
  "emens_indices"
)

## ---- dmsfsc-comparison-sensitivities
comparison_table(
  bind_rows(lapply(
    c("clock", "phase"),
    function(axis) read_scenarios(paste0("tables/shapley_", axis, ".csv"))
  )),
  "shapley"
)
comparison_table(
  read_scenarios("tables/sensitivity_summary.csv"),
  "model_sensitivities"
)
comparison_table(read_scenarios("tables/meq_model_fit.csv"), "meq_fit")
comparison_table(
  read_scenarios("tables/meq_candidate_windows.csv"),
  "meq_windows"
)
chest <- bind_rows(lapply(
  c("clock", "phase"),
  function(axis)
    read_scenarios(paste0("tables/chest_", axis, "_significant_windows.csv"))
))
comparison_table(chest, "chest_windows")
comparison_table(
  bind_rows(lapply(
    c("clock", "phase"),
    function(axis)
      read_scenarios(paste0(
        "tables/primary_",
        axis,
        "_individual_timing.csv"
      )) |>
        mutate(axis = axis)
  )),
  "individual_timing"
)
comparison_table(
  read_scenarios("tables/matched_placement_summary.csv"),
  "matched_placement"
)
noncyclic <- bind_rows(lapply(
  names(paths)[-1],
  function(scenario)
    bind_rows(lapply(
      c("clock", "phase"),
      function(axis)
        read.csv(file.path(
          paths[[scenario]],
          paste0(
            "results/tables/noncyclic_sz_",
            axis,
            "_significant_windows.csv"
          )
        )) |>
          mutate(scenario = scenario)
    ))
))
comparison_table(noncyclic, "noncyclic_windows")

## ---- dmsfsc-comparison-figures
figure_contrast <- contrasts |>
  filter(contrast == "evening_morning_Extreme_Overall") |>
  mutate(
    axis = factor(
      axis,
      levels = c("clock", "phase"),
      labels = c("Local clock time", "Relative midsleep time")
    ),
    scenario = factor(scenario, levels = names(paths))
  )
write.csv(
  figure_contrast,
  file.path(root, "data/intermediate/figure_contrasts_data.csv"),
  row.names = FALSE
)
p <- ggplot(
  figure_contrast,
  aes(x, ratio, colour = scenario, fill = scenario)
) +
  geom_hline(yintercept = 1, colour = "grey55") +
  geom_ribbon(
    aes(ymin = ratio_lower, ymax = ratio_upper),
    alpha = .08,
    colour = NA
  ) +
  geom_line(linewidth = .7) +
  facet_wrap(~axis, scales = "free_x", ncol = 1) +
  scale_y_log10(breaks = c(.01, .1, 1, 10, 100)) +
  scale_colour_manual(values = c("grey45", "#6F42A6", "#00856A"), labels = scenario_labels) +
  scale_fill_manual(values = c("grey45", "#6F42A6", "#00856A"), labels = scenario_labels) +
  labs(
    x = "Time (h)",
    y = "Extreme evening / morning exposure ratio",
    colour = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(
  file.path(root, "comparison/contrasts.png"),
  p,
  width = 10,
  height = 9,
  dpi = 180,
  bg = "white"
)
ggsave(file.path(root, "comparison/contrasts.pdf"), p, width = 10, height = 9)
figure_daylength <- bind_rows(lapply(
  names(paths),
  function(scenario)
    mutate(
      read.csv(file.path(
        paths[[scenario]],
        "results/tables/daylight_contrasts.csv"
      )),
      scenario = scenario
    )
)) |>
  filter(contrast == "Longer_minus_shorter") |>
  mutate(scenario = factor(scenario, levels = names(paths)))
write.csv(
  figure_daylength,
  file.path(root, "data/intermediate/figure_daylength_data.csv"),
  row.names = FALSE
)
p <- ggplot(
  figure_daylength,
  aes(x, ratio, colour = scenario, fill = scenario)
) +
  geom_hline(yintercept = 1, colour = "grey55") +
  geom_ribbon(
    aes(ymin = pointwise_ratio_lower, ymax = pointwise_ratio_upper),
    alpha = .08,
    colour = NA
  ) +
  geom_line(linewidth = .7) +
  scale_y_log10() +
  scale_colour_manual(values = c("grey45", "#6F42A6", "#00856A"), labels = scenario_labels) +
  scale_fill_manual(values = c("grey45", "#6F42A6", "#00856A"), labels = scenario_labels) +
  labs(
    x = "Time relative to the scenario midsleep metric (h)",
    y = "Longer / shorter ratio of ratios",
    colour = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(
  file.path(root, "comparison/daylength.png"),
  p,
  width = 10,
  height = 6,
  dpi = 180,
  bg = "white"
)
ggsave(file.path(root, "comparison/daylength.pdf"), p, width = 10, height = 6)
figure_prc <- prc |>
  filter(context == "Overall") |>
  mutate(
    prc_baseline = factor(
      prc_baseline,
      levels = c("Published raw", "Drift corrected")
    ),
    scenario = factor(scenario, levels = rev(names(paths)))
  )
write.csv(
  figure_prc,
  file.path(root, "data/intermediate/figure_prc_data.csv"),
  row.names = FALSE
)
p <- ggplot(figure_prc, aes(estimate, scenario)) +
  scale_y_discrete(labels = scenario_labels) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = .12,
    orientation = "y"
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~prc_baseline) +
  labs(x = "Extreme evening minus morning PRC-weighted index", y = NULL) +
  theme_bw(base_size = 11)
ggsave(
  file.path(root, "comparison/prc.png"),
  p,
  width = 10,
  height = 4,
  dpi = 180,
  bg = "white"
)
ggsave(file.path(root, "comparison/prc.pdf"), p, width = 10, height = 4)

## ---- dmsfsc-verification
checks <- list()
add_check <- function(check, passed) {
  checks[[length(checks) + 1L]] <<- tibble(check, passed = isTRUE(passed))
}
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
  add_check(
    paste(
      placement,
      "matched observation rows, response, day labels and AR ordering"
    ),
    identical(
      matched$frames[[placement]][unchanged],
      changed$frames[[placement]][unchanged]
    )
  )
}
folds <- lapply(
  paths[-1],
  function(path)
    read.csv(file.path(path, "results/validation/participant_folds.csv"))
)
add_check(
  "Identical held-out participant folds",
  identical(folds[[1]], folds[[2]])
)
for (scenario in names(paths)[-1]) {
  prep <- readRDS(file.path(paths[[scenario]], "results/data/prepared.rds"))
  registry <- read.csv(file.path(
    paths[[scenario]],
    "results/validation/participant_folds.csv"
  ))
  for (i in 1:5) {
    training <- filter(
      prep$frames$glasses,
      !participant_key %in% registry$participant_key[registry$fold == i]
    )
    expected <- fit_preprocess(training)
    for (axis in c("clock", "phase")) {
      out <- readRDS(file.path(
        paths[[scenario]],
        paste0("results/validation/new_participant_fold_", i, "_", axis, ".rds")
      ))
      add_check(
        paste(
          scenario,
          axis,
          "fold",
          i,
          "training-only cutoffs and photoperiod centers"
        ),
        isTRUE(all.equal(out$cutoffs, expected$cutoffs)) &&
          isTRUE(all.equal(out$centers, expected$centers))
      )
    }
  }
}
stage_checks <- bind_rows(lapply(
  names(paths)[-1],
  function(scenario)
    read.csv(file.path(
      paths[[scenario]],
      "results/provenance/scenario_stage_checks.csv"
    )) |>
      mutate(scenario = scenario)
))
comparison_table(stage_checks, "scenario_stage_checks")
add_check(
  "All scenario-stage completion, Shapley and noncyclic-basis checks pass",
  all(stage_checks$passed)
)
stages <- c(
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
)
stage_runs <- bind_rows(lapply(names(paths)[-1], function(scenario) {
  bind_rows(lapply(stages, function(stage) {
    read.csv(file.path(
      paths[[scenario]],
      "results/provenance",
      paste0(stage, "_run.csv")
    ))
  })) |>
    mutate(scenario = scenario)
}))
comparison_table(stage_runs, "completed_stages")
source_files <- list.files(
  file.path(root, "source-snapshot"),
  full.names = TRUE,
  pattern = "[.]R$"
)
for (scenario in names(paths)[-1]) {
  scenario_sources <- file.path(paths[[scenario]], "R", basename(source_files))
  add_check(
    paste(
      scenario,
      "frozen scientific producers match the initial source snapshot"
    ),
    identical(
      unname(tools::md5sum(source_files)),
      unname(tools::md5sum(scenario_sources))
    )
  )
}
original_hashes <- read.csv(file.path(root, "original-input-hashes.csv"))
original_hashes$after_md5 <- unname(tools::md5sum(original_hashes$path))
original_hashes$unchanged <- original_hashes$md5 == original_hashes$after_md5
original_hashes$authorized_display_change <- FALSE
diary_source <- which(endsWith(
  original_hashes$path,
  "/R/sleep_diary_alignment.R"
))
without_plot <- function(path) {
  lines <- readLines(path)
  markers <- grep("^## ---- ", lines)
  starts <- grep("^## ---- sleep-diary-(plot|scatter-data)$", lines)
  stopifnot(length(starts) == 2L)
  omit <- unlist(lapply(starts, function(start) {
    following <- markers[markers > start]
    end <- if (length(following)) following[1] - 1L else length(lines)
    seq.int(start, end)
  }))
  lines[-omit]
}
display_only <- identical(
  without_plot(file.path(root, "source-snapshot/sleep_diary_alignment.R")),
  without_plot("R/sleep_diary_alignment.R")
)
original_hashes$authorized_display_change[diary_source] <-
  !original_hashes$unchanged[diary_source] & display_only
comparison_table(original_hashes, "original_input_preservation")
add_check(
  "Original scientific inputs and models unchanged; authorized diary display edits recorded",
  all(original_hashes$unchanged | original_hashes$authorized_display_change)
)
checks <- bind_rows(checks)
comparison_table(checks, "implementation_checks")
capture.output(sessionInfo(), file = file.path(root, "comparison-session.txt"))
stopifnot(all(checks$passed))
