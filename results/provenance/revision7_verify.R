# Verify the new inferential family and its propagation without refitting models.
source('R/setup.R')
source('R/models.R');source('R/inference.R');source('R/meq.R')
checks <- list()
check <- function(name,ok) {
  stopifnot(isTRUE(ok))
  checks[[length(checks)+1L]] <<- tibble(check=name,passed=TRUE)
}
read_table <- function(name) read.csv(paste0('results/tables/',name,'.csv'))
protected <- read.csv('results/provenance/revision7_fitted_inputs_before.csv')
check('All fitted models, prepared observations and validation inputs remain byte-identical',
  identical(unname(tools::md5sum(protected$file)),protected$md5))
for(label in c('primary','chest','noncyclic_sz'))for(axis in c('clock','phase')) {
  m <- read_table(paste0(label,'_',axis,'_inference_metadata'))
  check(paste(label,axis,'95% simultaneous within-coordinate family'),
    m$axis_family_coverage==.95 && m$draws==10000 && m$minimum_cell_support==4)
}
for(label in c('meq_defined','meq_matched_msfsc'))for(axis in c('clock','phase')) {
  m <- read_table(paste0(label,'_',axis,'_inference_metadata'))
  check(paste(label,axis,'95% within-coordinate classification sensitivity'),
    m$within_fit_coverage==.95 && m$draws==10000 && m$minimum_cell_support==4)
}
b <- readRDS('results/models/primary_phase.rds')
inf <- readRDS('results/data/primary_phase_inference.rds')
av <- read_table('primary_phase_interval_averages')
V <- covariance_info(b$fit)$V
for(i in seq_len(nrow(av))) {
  ix <- which(wrap24(inf$cell_design$x-av$start[i])<av$duration_hours[i]-1e-10)
  L <- colMeans((inf$pattern_matrices$Q4_Overall-inf$pattern_matrices$Q1_Overall)[ix,,drop=FALSE])
  eta <- sum(L*coef(b$fit));se <- sqrt(drop(L%*%V%*%L))
  expected <- 10^(eta+c(0,-1,1)*qnorm(.975)*se)
  check(paste('Primary window',i,'full-covariance average and interval'),
    max(abs(expected-as.numeric(av[i,c('ratio','ratio_lower','ratio_upper')])))<1e-10)
}
windows <- inf$windows |> filter(contrast=='evening_morning_Extreme_Overall') |> arrange(start)
check('Two extreme-group phase intervals lie before the DLMO guide and after MSFsc',
  nrow(windows)==2 && windows$end[1]< -6.95 && windows$start[2]>0)
for(prefix in c('daylight','daylight_quartiles','daylight_p10_p90')) {
  z <- read_table(paste0(prefix,'_primary_interval_summaries'))
  w <- z |> distinct(start_phase,end_phase) |> arrange(start_phase)
  check(paste(prefix,'summaries use the updated primary intervals'),
    max(abs(w$start_phase-windows$start),abs(w$end_phase-windows$end))<1e-12)
  point <- read_table(paste0(prefix,'_pointwise_contrasts'))
  check(paste(prefix,'95% pointwise bands retain their original definition'),
    max(abs(point$pointwise_lower-(point$estimate-qnorm(.975)*point$se)),
        abs(point$pointwise_upper-(point$estimate+qnorm(.975)*point$se)))<1e-12)
}
civil <- read_table('daylight_civil_primary_intervals') |> distinct(start_phase,end_phase) |> arrange(start_phase)
check('Civil daylight summaries use the updated primary intervals',
  max(abs(civil$start_phase-windows$start),abs(civil$end_phase-windows$end))<1e-12)
figwindows <- read.csv('data/intermediate/manuscript/primary-phase-intervals.csv')
check('Manuscript areas of interest match the primary inference',
  max(abs(figwindows$start-windows$start),abs(figwindows$end-windows$end))<1e-12)
fig <- read.csv('data/intermediate/manuscript/figure-2-patterns-contrasts.csv')
source <- bind_rows(lapply(c('clock','phase'),function(axis)
  read_table(paste0('primary_',axis,'_contrasts')) |> mutate(axis=axis))) |>
  filter(context=='Overall',family=='Evening minus morning') |>
  arrange(axis,contrast,x)
fig <- fig |> arrange(axis,contrast,x)
check('Figure 2 intervals and exclusions match the current inference',
  nrow(fig)==nrow(source) && max(abs(fig$ratio_lower-source$ratio_lower),
    abs(fig$ratio_upper-source$ratio_upper))<1e-10 && identical(fig$significant,source$significant))
rm(b,inf);gc()
prepared <- readRDS('results/data/prepared.rds')
raw <- meq_prepare(prepared)
meq_verify(raw,read_table('meq_all_curves'),read_table('meq_all_contrasts'),prepared$cutoffs)
save_table(bind_rows(checks),'revision7_verification','results/provenance')
capture.output(sessionInfo(),file='results/provenance/revision7_verification_runtime.txt')
print(bind_rows(checks),n=Inf,width=Inf)
