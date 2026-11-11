# Standalone data imports; prefer a locally restored package library.
project_root <- normalizePath(getwd(), mustWork = TRUE)
upstream_root <- normalizePath(Sys.getenv('MELIDOS_MAIN', file.path(project_root, '..', 'ZaunerEtAl_bioRxiv_2026')), mustWork = FALSE)
local_library <- file.path(project_root, 'renv/library/standalone')
runtime_series <- paste0('R-',R.version$major,'.',strsplit(R.version$minor,'.',fixed=TRUE)[[1]][1])
existing_library <- Sys.glob(file.path(upstream_root, 'renv/library', '*', runtime_series, R.version$platform))
.libPaths(c(local_library, Sys.getenv('MELIDOS_R_LIBRARY', ''), existing_library, .libPaths()))
stopifnot(as.character(getRversion()) == '4.6.1')
omp <- Sys.getenv('MELIDOS_OPENMP_LIBRARY', '/opt/homebrew/opt/libomp/lib/libomp.dylib')
if (file.exists(omp) && !any(grepl('libomp', names(getLoadedDLLs())))) dyn.load(omp, local = FALSE)
suppressPackageStartupMessages({library(mgcv); library(dplyr); library(tidyr); library(ggplot2)})
invisible(loadNamespace('gratia'))
options(dplyr.summarise.inform = FALSE, stringsAsFactors = FALSE)
for (d in c('results/data','results/models','results/subsets','results/figures','results/tables','results/diagnostics','results/validation','results/provenance','results/logs')) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cfg <- list(seed = 22102026L, k = 12L, kp = 10L, kphoto = 5L,
            draws = 10000L, bootstrap = 2000L, nthreads = 2L,
            support_n = 4L, step = 1/6, curve_coverage = .95, version = '1.0.0')
set.seed(cfg$seed)
count_sample <- function(d) tibble(participants = n_distinct(d$participant_key), participant_days = n_distinct(d$participant_day_key), half_hours = nrow(d))
sample_label <- function(d) { n <- count_sample(d); sprintf('%s participants; %s participant-days; %s half-hours', n$participants, n$participant_days, format(n$half_hours, big.mark = ',')) }
save_table <- function(x, name, dir = 'results/tables') { write.csv(x, file.path(dir, paste0(name,'.csv')), row.names = FALSE, na = ''); invisible(x) }
save_plot <- function(p, name, width = 11, height = 8) {
  ggsave(file.path('results/figures', paste0(name,'.png')), p, width = width, height = height, dpi = 180, bg = 'white')
  ggsave(file.path('results/figures', paste0(name,'.pdf')), p, width = width, height = height, bg = 'white')
  invisible(p)
}
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.spacing.x = grid::unit(16, 'pt'), plot.caption = element_text(hjust = 0), legend.position = 'bottom'))
clock_label <- function(x) { m <- round(x * 60); sprintf('%02d:%02d', (m %/% 60) %% 24, m %% 60) }
wrap24 <- function(x) x %% 24
wrap12 <- function(x) ((x + 12) %% 24) - 12
# Display transformations leave model-scale estimates and inference intact.
backtransform_intensity <- function(eta) 10^eta - 0.1
intensity_display <- function(z) {
  for(v in c('estimate','lower','upper','group_estimate')) if(v %in% names(z))
    z[[paste0(v,'_melEDI')]] <- backtransform_intensity(z[[v]])
  z
}
contrast_display <- function(z) {
  z$ratio <- 10^z$estimate
  z$ratio_lower <- 10^z$lower
  z$ratio_upper <- 10^z$upper
  z
}
window_display <- function(z) {
  z$minimum_ratio <- 10^z$minimum_difference
  z$maximum_ratio <- 10^z$maximum_difference
  z
}
atomic_rds <- function(x, path) { tmp <- paste0(path, '.', Sys.getpid(), '.tmp'); saveRDS(x, tmp, compress = FALSE); stopifnot(file.rename(tmp, path)); invisible(x) }
log_step <- function(...) { cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ..., '\n'); flush.console() }
