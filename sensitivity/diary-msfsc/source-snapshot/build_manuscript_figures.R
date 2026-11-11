# Compose publication figures from the saved chronotype analysis.

#####
# Step 1: Load saved results and display helpers
#####

source("R/setup.R")
source("R/plots.R")
library(patchwork)

protected_paths <- c(list.files('results/models', full.names = TRUE),
  list.files('results/data', full.names = TRUE), 'analysis.qmd', 'references.bib')
protected_paths <- protected_paths[!dir.exists(protected_paths)]
protected_before <- tibble(path = protected_paths,
  md5 = unname(tools::md5sum(protected_paths)))

figure_dir <- "results/figures/manuscript"
data_dir <- "data/intermediate/manuscript"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- character()
read_result <- function(name) {
  path <- file.path("results/tables", paste0(name, ".csv"))
  inputs <<- c(inputs, path)
  tibble::as_tibble(read.csv(path))
}
export_data <- function(data, name) {
  write.csv(data, file.path(data_dir, paste0(name, ".csv")), row.names = FALSE)
}
save_pair <- function(plot, stem, width, height, title, caption) {
  for (ext in c("png", "pdf")) {
    ggsave(file.path(figure_dir, paste0(stem, ".", ext)), plot,
      width = width, height = height, dpi = 300, bg = "white")
  }
  wrapped <- paste(strwrap(caption, width = 115), collapse = "\n")
  standalone <- plot + plot_annotation(title = title, caption = wrapped,
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
      plot.caption = element_text(size = 8.5, hjust = 0, lineheight = 1.1)))
  extra_height <- 0.4 + 0.15 * length(strsplit(wrapped, "\n")[[1]])
  for (ext in c("png", "pdf")) {
    ggsave(file.path(figure_dir, paste0(stem, "-standalone.", ext)), standalone,
      width = width, height = height + extra_height, dpi = 300, bg = "white")
  }
}

publication_theme <- theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
  panel.spacing = grid::unit(20, "pt"), strip.background = element_rect(fill = "grey95"),
  strip.text = element_text(size = 10), legend.position = "bottom",
  legend.title = element_blank(), legend.text = element_text(size = 10),
  axis.title = element_text(size = 10), axis.text = element_text(size = 9),
  plot.tag = element_text(face = "bold", size = 13),
  plot.margin = margin(5, 5, 5, 5))
theme_set(publication_theme)
group_colors <- c(Morning = "#0072B2", Evening = "#D55E00")
extreme_colors <- c(Q1 = "#0072B2", Q4 = "#D55E00")
extreme_labels <- c(Q1 = "Extreme morning", Q4 = "Extreme evening")
phase_axis <- function(breaks = seq(-12, 12, 6), expand_mult = 0.025) {
  scale_x_continuous(breaks = breaks, limits = c(-12, 12),
    expand = expansion(mult = expand_mult))
}
windows <- read_result("primary_phase_significant_windows") |>
  filter(contrast == "evening_morning_Extreme_Overall") |> arrange(start)
export_data(windows, "primary-phase-intervals")
export_data(tibble(start = -6.95, end = -5.95, central = -6.45), "dlmo-guide")

#####
# Step 2: Chronotype distributions
#####

people <- read_result("chronotype_distribution_registry")
groups <- read_result("chronotype_distribution_msfsc_groups")
meq_groups <- read_result("chronotype_distribution_meq_groups")
cutoffs <- read_result("quartile_cutoffs")
people <- people |> mutate(quartile = factor(quartile, levels = paste0("Q", 1:4)))
quartile_colors <- c(Q1 = "#0072B2", Q2 = "#56B4E9", Q3 = "#E69F00", Q4 = "#D55E00")
quartile_labels <- setNames(paste0(groups$group, "\nn = ", groups$participants,
  "; ", groups$participant_days, " days"), groups$quartile)
hist_msf <- ggplot(people, aes(msf_hours, fill = quartile)) +
  geom_histogram(binwidth = 20 / 60, boundary = 0, closed = "left",
    color = "white", linewidth = 0.2) +
  geom_vline(xintercept = cutoffs$msf_hours, linetype = "dashed",
    color = "grey35", linewidth = 0.35) +
  scale_fill_manual(values = quartile_colors, labels = quartile_labels) +
  scale_x_continuous(breaks = 1:8, labels = clock_label) +
  scale_y_continuous(breaks = scales::breaks_width(5), expand = expansion(mult = c(0, 0.06))) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(tag = "A", x = "MSFsc (local clock time)", y = "Participants", fill = NULL)
meq_labels <- meq_groups |> mutate(
  x = c(23.25, 36, 50, 64, 77.75),
  label = paste0(sub(" ", "\n", sub(" type$", "", meq_type)), "\nn = ", participants))
meq_height <- max(table(people$meq))
hist_meq <- ggplot(people, aes(meq, fill = quartile)) +
  geom_histogram(binwidth = 1, boundary = 0.5, closed = "left",
    color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(30.5, 41.5, 58.5, 69.5), linetype = "dashed",
    color = "grey35", linewidth = 0.35) +
  geom_text(data = meq_labels, aes(x, y = meq_height * 1.42, label = label),
    inherit.aes = FALSE, size = 2.65, vjust = 1, lineheight = 0.95) +
  scale_fill_manual(values = quartile_colors, labels = quartile_labels) +
  scale_x_continuous(breaks = c(16, 30, 41, 58, 69, 86)) +
  coord_cartesian(xlim = c(16, 86)) +
  scale_y_continuous(breaks = scales::breaks_width(4), expand = expansion(mult = c(0, 0.04))) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(tag = "B", x = "MEQ score", y = "Participants", fill = NULL)
figure1 <- wrap_plots(hist_msf, hist_meq, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")
export_data(people |> select(msf_hours, quartile, meq, meq_type), "figure-1-chronotypes-data")
export_data(groups, "figure-1-chronotypes-groups")
export_data(meq_labels, "figure-1-chronotypes-meq-labels")
export_data(cutoffs, "figure-1-chronotypes-cutoffs")
msf_bins <- ggplot2::layer_data(hist_msf, 1)
meq_bins <- ggplot2::layer_data(hist_meq, 1)
export_data(msf_bins, "figure-1-chronotypes-msfsc-bins")
export_data(meq_bins, "figure-1-chronotypes-meq-bins")
stopifnot(sum(msf_bins$count) == 140, sum(meq_bins$count) == 140)
save_pair(figure1, "figure-1-chronotypes", 6.6, 6.7,
  "Chronotype distribution in the MeLiDos sample",
  "A: corrected midsleep on free days (MSFsc), in 20-minute bins; dashed lines mark the quartiles. B: MEQ scores, in one-point bins, with category boundaries and counts. Colors identify the same four MSFsc quartiles in both panels. Each of the 140 participants contributes once. The shared legend gives participants and available participant-days.")

#####
# Step 3: Exposure patterns and aligned contrasts
#####

curves <- bind_rows(read_result("primary_clock_curves") |> mutate(axis = "clock"),
  read_result("primary_phase_curves") |> mutate(axis = "phase")) |>
  filter(context == "Overall") |> mutate(extremity = factor(extremity, c("Moderate", "Extreme")),
    direction = factor(direction, c("Morning", "Evening")))
contrasts <- bind_rows(read_result("primary_clock_contrasts") |> mutate(axis = "clock"),
  read_result("primary_phase_contrasts") |> mutate(axis = "phase")) |>
  filter(context == "Overall", family == "Evening minus morning") |>
  mutate(extremity = factor(sub("Evening minus morning: ", "", description), c("Moderate", "Extreme")))
export_data(curves, "figure-2-patterns-curves")
export_data(contrasts, "figure-2-patterns-contrasts")
panels <- list()
for (coordinate in c("clock", "phase")) {
  z <- filter(curves, axis == coordinate)
  d <- filter(contrasts, axis == coordinate) |>
    group_by(contrast) |> arrange(x, .by_group = TRUE) |>
    mutate(segment = cumsum(c(TRUE, diff(significant) != 0))) |> ungroup()
  horizontal_scale <- if (coordinate == "clock") {
    scale_x_continuous(breaks = seq(0, 24, 6), labels = function(v) ifelse(v == 24, "24:00", clock_label(v)),
      limits = c(0, 24), expand = expansion(mult = 0.025))
  } else phase_axis(expand_mult = 0.05)
  pattern <- ggplot(z, aes(x, estimate_melEDI, color = direction, fill = direction))
  contrast <- ggplot(d, aes(x, ratio))
  if (coordinate == "phase") {
    pattern <- pattern + dlmo_reference_band() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35)
    contrast <- contrast + dlmo_reference_band() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35)
  }
  pattern <- pattern +
    geom_ribbon(aes(ymin = lower_melEDI, ymax = upper_melEDI), alpha = 0.16, color = NA) +
    geom_line(linewidth = 0.6) + facet_wrap(~extremity, nrow = 1) +
    scale_color_manual(values = group_colors) + scale_fill_manual(values = group_colors) +
    horizontal_scale + intensity_scale(limits = range(c(0, curves$lower_melEDI, curves$upper_melEDI)),
      expand = expansion(mult = c(0.03, 0.15))) +
    labs(tag = if (coordinate == "clock") "A" else "B", x = NULL,
      y = "Melanopic EDI (lx)", color = NULL, fill = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (coordinate == "clock") {
    counts <- distinct(z, extremity, direction, participants) |>
      mutate(label = paste0(direction, ": n = ", participants), v = if_else(direction == "Morning", 1.1, 2.5))
    pattern <- pattern + geom_text(data = counts,
      aes(x = Inf, y = Inf, label = label, color = direction, vjust = v),
      inherit.aes = FALSE, hjust = 1.05, size = 2.8, show.legend = FALSE)
    export_data(counts, "figure-2-patterns-counts")
  }
  contrast <- contrast + geom_hline(yintercept = 1, color = "grey50", linewidth = 0.35) +
    geom_ribbon(aes(ymin = ratio_lower, ymax = ratio_upper), fill = "#6A3D9A", alpha = 0.18) +
    geom_line(color = "#6A3D9A", linewidth = 0.5) +
    geom_line(data = filter(d, significant), aes(group = interaction(contrast, segment)),
      color = "#6A3D9A", linewidth = 1.2) +
    facet_wrap(~extremity, nrow = 1) + horizontal_scale +
    scale_y_log10(limits = range(c(1, contrasts$ratio_lower, contrasts$ratio_upper)),
      breaks = c(0.1, 1, 10, 100), labels = as.character, oob = scales::oob_keep,
      expand = expansion(mult = c(0.03, 0.03))) +
    labs(x = if (coordinate == "clock") "Local clock time" else "Time relative to MSFsc (h)",
      y = "Evening / morning") + theme(strip.text = element_blank(), strip.background = element_blank())
  if (coordinate == "phase") {
    # Centre labels over the significant grid points. Place them in the clear
    # space above the pre-MSFsc band and below the post-MSFsc band.
    callouts <- windows |> mutate(label = if_else(end < 0, "pre-MSFsc", "post-MSFsc"),
      target_x = (start + end - cfg$step) / 2,
      extremity = factor("Extreme", levels = c("Moderate", "Extreme")))
    callouts$target_y <- vapply(callouts$target_x, function(v) {
      zz <- filter(d, extremity == "Extreme")
      zz$ratio[which.min(abs(zz$x - v))]
    }, numeric(1))
    callouts <- callouts |> mutate(label_x = target_x,
      label_y = if_else(end < 0, 5.0, 0.12),
      # Gaps are set on the plotted log scale at both ends of each leader.
      leader_start = label_y * 10^if_else(end < 0, -0.28, 0.28),
      leader_end = target_y * 10^if_else(end < 0, 0.18, -0.18))
    contrast <- contrast +
      geom_segment(data = callouts,
        aes(x = label_x, xend = target_x, y = leader_start, yend = leader_end),
        inherit.aes = FALSE, linewidth = 0.3, color = "black") +
      geom_text(data = callouts, aes(x = label_x, y = label_y, label = label),
        inherit.aes = FALSE, size = 2.6, color = "black")
    export_data(callouts, "figure-2-patterns-phase-callouts")
  }
  panels <- c(panels, list(pattern, contrast))
}
figure2 <- wrap_plots(panels, ncol = 1, heights = c(2, 1, 2, 1), guides = "collect") &
  theme(legend.position = "bottom")
save_pair(figure2, "figure-2-patterns", 7.2, 8.6,
  "Chronotype differences persist after alignment to corrected midsleep",
  "A: local clock time. B: time relative to MSFsc. Each block pairs exposure curves with shorter evening/morning contrast panels. Bands are 95% simultaneous within each coordinate: across the 12 group-by-day-context exposure curves above, and the 20 registered contrasts below. One denotes no group difference; bold lavender segments mark its supported exclusion. Black labels and vertical leaders identify the pre-MSFsc and post-MSFsc intervals. Dashed vertical lines mark MSFsc at zero. Ratios refer to melanopic EDI + 0.1. Stripes mark the external DLMO guide. Common photoperiod and work/free-day distributions; 140 participants and 778 days.")

#####
# Step 4: Daylength and civil daylight
#####

day_curves <- read_result("daylight_patterns") |>
  mutate(scenario = factor(scenario, c("Shorter", "Typical", "Longer")))
day_contrasts <- read_result("daylight_pointwise_contrasts") |>
  group_by(contrast) |> arrange(x, .by_group = TRUE) |>
  mutate(segment = cumsum(c(TRUE, diff(pointwise_significant) != 0))) |> ungroup()
day_ratios <- day_contrasts |> filter(contrast != "Longer_minus_shorter") |>
  mutate(scenario = factor(contrast, c("Shorter", "Typical", "Longer")))
day_change <- day_contrasts |> filter(contrast == "Longer_minus_shorter")
civil <- read_result("daylight_civil_availability") |>
  mutate(daylength = factor(daylength, c("Shorter", "Longer")))
export_data(day_curves, "figure-3-daylength-curves")
export_data(day_ratios, "figure-3-daylength-ratios")
export_data(day_change, "figure-3-daylength-moderation")
export_data(civil, "figure-3-daylength-civil")
targets <- read_result("daylight_scenario_targets")
export_data(targets, "figure-3-daylength-targets")
stopifnot(identical(as.numeric(targets$photo), c(-1, 0, 1)))
scenario_labels <- c(Shorter = "Shorter days (-1 h)", Typical = "Typical days (0 h)", Longer = "Longer days (+1 h)")
day_a <- ggplot(day_curves, aes(x, estimate_melEDI, color = quartile, fill = quartile)) +
  primary_phase_interval_band(boundaries = FALSE) + dlmo_reference_band() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lower_melEDI, ymax = upper_melEDI), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.6) + facet_wrap(~scenario, ncol = 3, labeller = as_labeller(scenario_labels)) +
  scale_color_manual(values = extreme_colors, labels = extreme_labels) +
  scale_fill_manual(values = extreme_colors, labels = extreme_labels) +
  scale_y_continuous(transform = LightLogR::symlog_trans(base = 10, thr = 1, scale = 1),
    breaks = c(0, 1, 10, 100, 1000)) + phase_axis() +
  labs(tag = "A", x = "Time relative to MSFsc (h)", y = "Melanopic EDI (lx)", color = NULL, fill = NULL)
# The same pointwise bands and support masks as in the analysis supplement.
ratio_panels <- list()
for (panel_data in list(day_ratios, day_change)) {
  limits <- range(c(1, panel_data$ratio, panel_data$pointwise_ratio_lower[panel_data$supported],
    panel_data$pointwise_ratio_upper[panel_data$supported]))
  p <- ggplot(panel_data, aes(x, ratio)) + primary_phase_interval_band(boundaries = FALSE) + dlmo_reference_band() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = 1, color = "grey50", linewidth = 0.35) +
    geom_ribbon(aes(ymin = ifelse(supported, pointwise_ratio_lower, NA),
      ymax = ifelse(supported, pointwise_ratio_upper, NA)), fill = "#6A3D9A", alpha = 0.18) +
    geom_line(color = "#6A3D9A", linewidth = 0.5) +
    geom_line(data = filter(panel_data, pointwise_significant),
      aes(group = interaction(contrast, segment)), color = "#6A3D9A", linewidth = 1.2) +
    geom_point(data = filter(panel_data, !supported), shape = 4, size = 0.8, color = "grey50") +
    scale_y_log10(limits = limits, oob = scales::oob_keep,
      breaks = c(0.03, 0.1, 0.3, 1, 3, 10, 30, 100), labels = function(v) paste0(v, "×")) +
    phase_axis() + labs(x = "Time relative to MSFsc (h)")
  ratio_panels <- c(ratio_panels, list(p))
}
day_b <- ratio_panels[[1]] + facet_wrap(~scenario, ncol = 3, labeller = as_labeller(scenario_labels)) +
  labs(tag = "B", y = "Evening / morning")
day_c <- ratio_panels[[2]] + labs(tag = "C", y = "Longer / shorter", x = "Time relative to\nMSFsc (h)")
day_d <- ggplot(civil, aes(x, availability, color = quartile)) +
  primary_phase_interval_band(boundaries = FALSE) + dlmo_reference_band() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_line(linewidth = 0.6, show.legend = FALSE) +
  facet_wrap(~daylength, labeller = as_labeller(c(Shorter = "Shorter absolute days", Longer = "Longer absolute days"))) +
  scale_color_manual(values = extreme_colors, labels = extreme_labels) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  phase_axis(c(-12, 0, 12)) + labs(tag = "D", x = "Time relative to MSFsc (h)", y = "Civil daylight") +
  theme(legend.position = "none", strip.text = element_text(size = 8.5))
figure3 <- wrap_plots(A = day_a, B = day_b, C = day_c, D = day_d,
  design = "AAA\nBBB\nCDD", heights = c(1, 1, 1.1), guides = "collect") &
  theme(legend.position = "bottom")
save_pair(figure3, "figure-3-daylength", 7.6, 8.2,
  "Daylength moderates the post-MSFsc chronotype contrast",
  "A: extreme-group exposure curves at -1, 0 and +1 h of site-centered photoperiod. B: evening/morning ratios. C: longer/shorter ratio of those ratios. Bands are pointwise 95% intervals; bold sections exclude one with local support. D: civil daylight availability, using an absolute-daylength split at 13.25 h and equal participant and within-participant day weights. Borderless grey windows identify the primary exposure intervals; dashed vertical lines mark MSFsc at zero and stripes mark the DLMO guide. Exposure model: 140 participants, 778 days. Civil daylight: 70 extreme-group participants, 390 days.")

#####
# Step 5: PRC source and weighted indices
#####

prc_path <- "results/data/phase_response_results.rds"
inputs <- c(inputs, prc_path)
prc_result <- readRDS(prc_path)
offset <- -6.45
phase <- seq(-12, 12, by = 1 / 60)
kernel <- tibble(phase_msfsc = phase, published = prc_result$prc$predict(phase - offset),
  drift_corrected = published + 0.54) |>
  pivot_longer(c(published, drift_corrected), names_to = "curve", values_to = "response")
points <- as_tibble(prc_result$prc$points) |>
  mutate(phase_msfsc = wrap12(hours_from_dlmo_midpoint + offset))
indices <- prc_result$summary |> filter(dlmo_offset == offset) |>
  mutate(prc_baseline = factor(prc_baseline, c("Published raw", "Drift corrected")),
    context = factor(context, c("Overall", "Work", "Free")), direction = factor(direction, c("Morning", "Evening")))
export_data(kernel, "figure-4-prc-kernel")
export_data(points, "figure-4-prc-source-points")
export_data(indices, "figure-4-prc-indices")
stopifnot(nrow(points) == 18, nrow(indices) == 12,
  all(indices$extremity == "Extreme"), all(indices$exposure_scale == "Zero-exposure centered"))
kernel_labels <- c(published = "Published raw", drift_corrected = "Drift corrected")
ylimits <- range(c(kernel$response, points$phase_shift_hours)) + c(-0.12, 0.65)
prc_window_labels <- windows |> transmute(
  x = (start + end) / 2, y = ylimits[2] - 0.12,
  label = if_else(end < 0, "pre-MSFsc", "post-MSFsc"))
export_data(prc_window_labels, "figure-4-prc-window-labels")
prc_a <- ggplot(kernel, aes(phase_msfsc, response, color = curve, linetype = curve)) +
  primary_phase_interval_band(boundaries = FALSE) + dlmo_reference_band() +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = offset, color = "#8B6423", linewidth = 0.45, linetype = 3) +
  geom_vline(xintercept = 0, color = "#0072B2", linewidth = 0.45, linetype = 3) +
  geom_point(data = points, aes(phase_msfsc, phase_shift_hours), inherit.aes = FALSE,
    color = "grey50", alpha = 0.75, size = 1.4) + geom_line(linewidth = 0.7) +
  annotate("label", x = offset, y = ylimits[2] - 0.12, label = "DLMO guide\n-6.45 h",
    color = "#8B6423", fill = "white", linewidth = 0, vjust = 1, size = 2.9) +
  annotate("label", x = 0, y = ylimits[2] - 0.12, label = "MSFsc\n0 h",
    color = "#0072B2", fill = "white", linewidth = 0, vjust = 1, size = 2.9) +
  geom_label(data = prc_window_labels, aes(x, y, label = label), inherit.aes = FALSE,
    color = "grey30", fill = "white", linewidth = 0, vjust = 1, size = 2.9,
    label.padding = grid::unit(0.07, "lines")) +
  scale_color_manual(values = c(drift_corrected = "#6A3D9A", published = "#666666"), breaks = names(kernel_labels), labels = kernel_labels) +
  scale_linetype_manual(values = c(drift_corrected = "solid", published = "longdash"), breaks = names(kernel_labels), labels = kernel_labels) +
  scale_x_continuous(breaks = seq(-12, 12, 3), limits = c(-12, 12),
    sec.axis = sec_axis(~. - offset, breaks = seq(-6, 18, 6), name = "Pulse midpoint relative to DLMO (h)")) +
  scale_y_continuous(limits = ylimits, expand = expansion(mult = 0.01)) +
  labs(tag = "A", x = "Pulse midpoint relative to MSFsc (h)", y = "Phase response (h)", color = NULL, linetype = NULL)
prc_b <- ggplot(indices, aes(context, estimate, color = direction, group = direction)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2, linewidth = 0.35) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, position = position_dodge(0.45)) +
  geom_point(position = position_dodge(0.45), size = 2.2) + facet_wrap(~prc_baseline, nrow = 1) +
  scale_color_manual(values = group_colors, labels = c(Morning = "Extreme morning", Evening = "Extreme evening")) +
  labs(tag = "B", x = NULL, y = "PRC-weighted index", color = NULL)
emens <- read_result("emens_balance_summary") |>
  filter(dlmo_offset == offset) |>
  mutate(context = factor(context, c("Overall", "Work", "Free")),
    direction = factor(direction, c("Morning", "Evening")),
    method = "Broad-region balance")
export_data(emens, "figure-4-emens-balances")
stopifnot(nrow(emens) == 6L, all(emens$quartile %in% c("Q1", "Q4")))
prc_c <- ggplot(emens, aes(context, estimate, color = direction, group = direction)) +
  geom_hline(yintercept = 1, color = "grey55", linetype = 2, linewidth = 0.35) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, position = position_dodge(0.45)) +
  geom_point(position = position_dodge(0.45), size = 2.2) + facet_wrap(~method, nrow = 1) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = group_colors, labels = c(Morning = "Extreme morning", Evening = "Extreme evening")) +
  labs(tag = "C", x = NULL, y = "Advance / delay exposure", color = NULL)
prc_bottom <- wrap_plots(prc_b, prc_c, nrow = 1, widths = c(2.05, 1), guides = "collect") &
  theme(legend.position = "bottom")
figure4 <- wrap_plots(prc_a, prc_bottom, ncol = 1, heights = c(1.3, 1))
save_pair(figure4, "figure-4-prc", 7.2, 7.4,
  "PRC weighting and advance-to-delay exposure balance",
  "A: one-hour bright-white-light PRC from St Hilaire et al. (2012), reconstructed from 18 coordinates deposited by Walch et al. (2025). Dashed: published raw baseline; solid: drift corrected by +0.54 h. DLMO is aligned 6.45 h before MSFsc. B: extreme-group indices with the published raw PRC on the left and drift-corrected PRC on the right, using log10(melanopic EDI + 0.1) + 1. C: Emens-style sensitivity, dividing mean log10(fitted melanopic EDI + 1) in the advance region (DLMO -16 to -4 h) by the delay-region mean (-4 to +8 h). The line at one in C marks equal regional means within a group. Bars in B and C are 95% coefficient intervals; C uses percentiles of 10,000 joint draws. Grey windows are labelled pre-MSFsc and post-MSFsc; stripes mark the DLMO guide. Extreme groups: 70 participants, 390 days. These balances are not phase shifts in hours.")

#####
# Step 6: Record inputs and verify preservation
#####

source_paths <- unique(c(inputs, "R/build_manuscript_figures.R", "R/setup.R", "R/plots.R"))
write.csv(tibble(path = source_paths, md5 = unname(tools::md5sum(source_paths))),
  "results/provenance/manuscript_figure_inputs.csv", row.names = FALSE)
protected_after <- protected_before |> mutate(md5 = unname(tools::md5sum(path)))
stopifnot(identical(protected_before, protected_after))
write.csv(protected_before, "results/provenance/manuscript_figure_protected_before.csv", row.names = FALSE)
write.csv(protected_after, "results/provenance/manuscript_figure_protected_after.csv", row.names = FALSE)
writeLines(c("Command: Rscript --vanilla R/build_manuscript_figures.R",
  paste("Executed:", Sys.time()), "All saved models, data caches, analysis.qmd and references.bib unchanged.",
  "Both histogram panels retain all 140 participants; source PRC retains 18 responses.",
  "Figure data are paired CSVs in data/intermediate/manuscript.", capture.output(sessionInfo())),
  "results/provenance/manuscript_figure_runtime.txt")
cat("Four manuscript figures and four standalone counterparts exported as PNG and PDF.\n")
