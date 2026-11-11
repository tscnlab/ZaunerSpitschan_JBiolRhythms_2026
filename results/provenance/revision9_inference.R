# Primary exposure-curve simultaneous bands, retaining the fitted models,
# registered contrast inference, and pointwise daylength displays.
# Command: Rscript --vanilla results/provenance/revision9_inference.R
source('R/setup.R')
for(f in c('data','models','inference','plots','phase_response'))source(paste0('R/',f,'.R'))
protected<-c(list.files('results/models',full.names=TRUE,pattern='\\.rds$'),
  list.files('results/data',full.names=TRUE,pattern='\\.rds$'),
  list.files('results/tables',full.names=TRUE,pattern='^(daylight|phase_response).*\\.csv$'),
  'results/validation/held_out_predictions.csv')
protected<-protected[!grepl('results/data/primary_(clock|phase)_inference.rds$',protected)]
before<-tools::md5sum(protected)
save_table(tibble(path=names(before),md5=unname(before)),
  'revision9_protected_before','results/provenance')
curves<-checks<-list()
for(axis in c('clock','phase')) {
  stem<-paste0('primary_',axis)
  bundle<-readRDS(paste0('results/models/',stem,'.rds'))
  old<-readRDS(paste0('results/data/',stem,'_inference.rds'))
  inf<-infer_curves(bundle)
  stopifnot(identical(inf$contrasts,old$contrasts),identical(inf$windows,old$windows),
    identical(inf$curves$estimate,old$curves$estimate),
    identical(inf$curves$se,old$curves$se),
    identical(inf$curves$support_n,old$curves$support_n),
    inf$pattern_critical>qnorm(.975),
    identical(inf$pattern_critical,unname(quantile(inf$pattern_simulated_maxima,.95,type=8))),
    all(inf$curves$lower<=old$curves$lower+1e-12),
    all(inf$curves$upper>=old$curves$upper-1e-12))
  atomic_rds(inf,paste0('results/data/',stem,'_inference.rds'))
  curves[[axis]]<-mutate(inf$curves,axis=axis)
  checks[[axis]]<-tibble(axis=axis,curves=n_distinct(paste(inf$curves$quartile,inf$curves$context)),
    grid_points=n_distinct(inf$curves$x),coverage=inf$pattern_coverage,
    critical=inf$pattern_critical,draws=cfg$draws,contrasts_unchanged=TRUE,
    windows_unchanged=TRUE,estimates_unchanged=TRUE,standard_errors_unchanged=TRUE)
  rm(bundle,old,inf);gc()
}
curves<-bind_rows(curves)
for(context in c('Overall','Work','Free'))save_plot(pattern_figure(curves,context,'Near-eye'),
  paste0('primary_patterns_',tolower(context)),12,11)
# Reorder the existing PRC display without recalculating any weighted index.
phase_response_figures(readRDS('results/data/phase_response_results.rds'))
save_table(bind_rows(checks),'revision9_exposure_band_checks','results/provenance')
scores<-read.csv('results/tables/validation_scores.csv')
save_table(scores,'revision9_verified_predictive_scores','results/provenance')
after<-tools::md5sum(protected)
stopifnot(identical(before,after))
save_table(tibble(path=names(before),before_md5=unname(before),after_md5=unname(after),
  unchanged=before==after),'revision9_preservation','results/provenance')
writeLines(c('Command: Rscript --vanilla results/provenance/revision9_inference.R',
  paste('Executed:',Sys.time()),
  'Inputs: saved primary clock and phase fits and inference; saved PRC results; validation_scores.csv.',
  '95% exposure-band family: all 12 quartile-by-context curves over 144 points within each coordinate.',
  'No model refitting, contrast or selected-window change, PRC recalculation, or daylength change.',
  capture.output(sessionInfo())),'results/provenance/revision9_runtime.txt')
print(bind_rows(checks),width=Inf)
print(scores[,c('target','axis','predictive_R2')],row.names=FALSE)
