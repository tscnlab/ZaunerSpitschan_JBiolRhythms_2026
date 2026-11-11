# Compare timing-only and timing-plus-classification changes on matched records.
# Run from the project root after the fixed-groups primary and daylight stages.
source("R/setup.R")
root <- "sensitivity/diary-msfsc"
for (module in c("data", "models", "inference")) {
  source(file.path(root, "source-snapshot", paste0(module, ".R")))
}
scenarios <- tibble(
  scenario = c("matched-msfsc", "dmsfsc-fixed-groups", "dmsfsc"),
  alignment = c("Questionnaire", "Diary", "Diary"),
  classification = c("Questionnaire", "Questionnaire", "Diary"),
  label = c(
    "Questionnaire timing + groups",
    "Diary timing + questionnaire groups",
    "Diary timing + groups"
  )
)
export_comparison <- function(x, name) {
  write.csv(
    x,
    file.path(root, "comparison", paste0(name, ".csv")),
    row.names = FALSE
  )
  invisible(x)
}
windows <- read.csv("results/tables/primary_phase_significant_windows.csv") |>
  filter(contrast == "evening_morning_Extreme_Overall")

## ---- fixed-groups-average-contrasts
intervals <- curves <- daylength <- fits <- diagnostics <- samples <- list()
model_data <- list()
for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[i]
  result_path <- file.path(root, scenario, "results")
  bundle <- readRDS(file.path(result_path, "models/primary_phase.rds"))
  inference <- readRDS(file.path(
    result_path,
    "data/primary_phase_inference.rds"
  ))
  covariance <- covariance_info(bundle$fit)
  model_data[[scenario]] <- bundle$data
  samples[[scenario]] <- bind_cols(scenarios[i, ], count_sample(bundle$data))
  fits[[scenario]] <- mutate(model_summary(bundle), scenario = scenario)
  for (context in c("Overall", "Work", "Free")) {
    id <- paste0("evening_morning_Extreme_", context)
    weights <- inference$registry$W[id, ]
    contrast_matrix <- Reduce(
      `+`,
      Map(function(x, w) x * w, inference$cell_design$X, weights)
    )
    for (j in seq_len(nrow(windows))) {
      inside <- wrap24(inference$cell_design$x - windows$start[j]) <
        windows$duration_hours[j] - 1e-10
      summary <- linear_summary(
        matrix(colMeans(contrast_matrix[inside, , drop = FALSE]), nrow = 1),
        coef(bundle$fit),
        covariance$V
      )
      support <- inference$contrasts |>
        filter(contrast == id, x %in% inference$cell_design$x[inside])
      intervals[[length(intervals) + 1L]] <- bind_cols(
        scenarios[i, ],
        summary
      ) |>
        mutate(
          context = context,
          phase_label = ifelse(windows$end[j] < 0, "Pre-MSFsc", "Post-MSFsc"),
          start = windows$start[j],
          end = windows$end[j],
          grid_points = sum(inside),
          ratio = 10^estimate,
          ratio_lower = 10^lower,
          ratio_upper = 10^upper,
          supported = all(support$supported),
          minimum_cell_participants = min(support$minimum_cell_participants)
        )
    }
  }
  curves[[scenario]] <- inference$contrasts |>
    filter(contrast == "evening_morning_Extreme_Overall") |>
    transmute(
      scenario = scenario,
      label = scenarios$label[i],
      x = wrap12(x),
      ratio,
      lower = ratio_lower,
      upper = ratio_upper,
      supported,
      significant,
      band = "95% simultaneous"
    )
  daylength[[scenario]] <- read.csv(file.path(
    result_path,
    "tables/daylight_interval_comparison.csv"
  )) |>
    filter(window_source == "Original primary fixed") |>
    mutate(
      scenario = scenario,
      alignment = scenarios$alignment[i],
      classification = scenarios$classification[i],
      label = scenarios$label[i]
    )
  diagnostics[[scenario]] <- bind_rows(lapply(
    c("primary", "daylight"),
    function(model) {
      read.csv(file.path(
        result_path,
        "diagnostics",
        paste0(model, "_phase_verdict.csv")
      )) |>
        mutate(scenario = scenario, model = model)
    }
  ))
  rm(bundle, inference, covariance)
  gc()
}
intervals <- export_comparison(
  bind_rows(intervals),
  "fixed_groups_exposure_intervals"
)
daylength <- export_comparison(
  bind_rows(daylength),
  "fixed_groups_daylength_intervals"
)
export_comparison(bind_rows(samples), "fixed_groups_samples")
export_comparison(bind_rows(fits), "fixed_groups_model_fit")
export_comparison(bind_rows(diagnostics), "fixed_groups_diagnostics")

## ---- fixed-groups-comparison-checks
questionnaire <- model_data[["matched-msfsc"]]
fixed <- model_data[["dmsfsc-fixed-groups"]]
diary <- model_data[["dmsfsc"]]
timing_fields <- c("msf_hours", "phase", "phase_display", "t")
group_fields <- c("quartile", "direction", "extremity", "group_label")
checks <- tibble(
  check = c(
    "Timing-only scenario preserves all non-timing fields",
    "Timing-only scenario uses exactly the diary alignment",
    "Final classification change preserves all non-group fields",
    "All fixed-window exposure comparisons meet the four-person support rule",
    "All main fixed-window daylength comparisons meet the four-person support rule"
  ),
  passed = c(
    identical(
      fixed[setdiff(names(fixed), timing_fields)],
      questionnaire[setdiff(names(questionnaire), timing_fields)]
    ),
    identical(fixed[timing_fields], diary[timing_fields]),
    identical(
      fixed[setdiff(names(fixed), group_fields)],
      diary[setdiff(names(diary), group_fields)]
    ),
    all(intervals$supported),
    all(daylength$supported[daylength$target_specification == "fixed"])
  )
)
exports <- file.path(root, "dmsfsc-fixed-groups/results/provenance")
input_hashes <- read.csv(file.path(exports, "input-hashes.csv"))
checks <- bind_rows(
  checks,
  tibble(
    check = "Fixed-group model inputs and producers have unchanged hashes",
    passed = identical(
      input_hashes$md5,
      unname(tools::md5sum(input_hashes$path))
    )
  )
)
export_comparison(checks, "fixed_groups_checks")
stopifnot(all(checks$passed))

## ---- fixed-groups-comparison-figure
moderation_curves <- bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
  read.csv(file.path(
    root,
    scenarios$scenario[i],
    "results/tables/daylight_contrasts.csv"
  )) |>
    filter(contrast == "Longer_minus_shorter") |>
    transmute(
      scenario = scenarios$scenario[i],
      label = scenarios$label[i],
      x = wrap12(x),
      ratio,
      lower = pointwise_ratio_lower,
      upper = pointwise_ratio_upper,
      supported,
      significant = pointwise_significant,
      band = "95% pointwise"
    )
}))
figure_data <- bind_rows(
  bind_rows(curves) |>
    mutate(comparison = "A  Extreme evening / morning exposure"),
  moderation_curves |>
    mutate(comparison = "B  Longer / shorter daylength moderation")
) |>
  mutate(label = factor(label, levels = scenarios$label)) |>
  arrange(comparison, label, x)
write.csv(
  figure_data,
  file.path(root, "data/intermediate/fixed_groups_comparison_data.csv"),
  row.names = FALSE
)
figure_limits <- figure_data |>
  summarise(ymin = min(lower) / 1.1, ymax = max(upper) * 1.1, .by = comparison)
window_background <- cross_join(windows, figure_limits)
plot <- ggplot(figure_data, aes(x, ratio, colour = label, fill = label)) +
  geom_rect(
    data = window_background,
    aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey75",
    alpha = .22,
    colour = NA
  ) +
  geom_hline(yintercept = 1, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .10, colour = NA) +
  geom_line(linewidth = .8) +
  facet_wrap(~comparison, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(-12, 12, 3)) +
  scale_y_log10(expand = expansion(mult = 0)) +
  scale_colour_manual(values = c("#6F42A6", "#B58900", "#00856A")) +
  scale_fill_manual(values = c("#6F42A6", "#B58900", "#00856A")) +
  labs(
    x = "Time relative to the scenario midsleep measure (h)",
    y = "Ratio",
    colour = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
ggsave(
  file.path(root, "comparison/fixed_groups.png"),
  plot,
  width = 9,
  height = 8,
  dpi = 180,
  bg = "white"
)
ggsave(
  file.path(root, "comparison/fixed_groups.pdf"),
  plot,
  width = 9,
  height = 8,
  bg = "white"
)
capture.output(
  sessionInfo(),
  file = file.path(root, "fixed-groups-comparison-session.txt")
)
print(
  daylength |>
    filter(
      target_specification == "fixed",
      contrast == "Longer_minus_shorter"
    ) |>
    select(
      label,
      phase_label,
      ratio,
      ratio_lower,
      ratio_upper,
      minimum_cell_participants
    )
)
