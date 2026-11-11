# Adopt 95% simultaneous inference within each time coordinate from saved fits.
# Command: Rscript --vanilla results/provenance/revision7_inference.R
source('R/setup.R')
for(f in c('data','models','inference','plots','diagnostics','validation','daylight','noncyclic_sz','meq'))
  source(paste0('R/',f,'.R'))
stopifnot(cfg$curve_coverage==.95,cfg$support_n==4L)
protected <- c(list.files('results/models',pattern='\\.rds$',full.names=TRUE),
  'results/data/prepared.rds','results/validation/participant_folds.csv',
  'results/validation/held_out_predictions.csv')
protected <- protected[!grepl('_inference\\.rds$',protected)]
before <- tools::md5sum(protected)
save_table(tibble(file=names(before),md5=unname(before)),
  'revision7_fitted_inputs_before','results/provenance')
old_tables <- list.files('results/tables',pattern='(significant_windows|candidate_windows|interval_averages|inference_metadata|primary_interval_summaries|civil_primary_intervals)\\.csv$',full.names=TRUE)
old_tables <- old_tables[!file.exists(file.path('results/provenance/revision7-before',basename(old_tables)))]
invisible(file.copy(old_tables,'results/provenance/revision7-before',overwrite=FALSE))
checks <- list()
for(label in c('primary','chest')) {
  allcontrasts <- list()
  for(axis in c('clock','phase')) {
    stem <- paste(label,axis,sep='_')
    b <- readRDS(paste0('results/models/',stem,'.rds'))
    old <- readRDS(paste0('results/data/',stem,'_inference.rds'))
    inf <- infer_curves(b)
    stopifnot(identical(inf$contrasts$estimate,old$contrasts$estimate),
      identical(inf$contrasts$se,old$contrasts$se),
      identical(inf$contrasts$supported,old$contrasts$supported),
      inf$critical<=old$critical+1e-12,
      all(!old$contrasts$significant|inf$contrasts$significant))
    checks[[stem]] <- tibble(model=stem,old_critical=old$critical,
      new_critical=inf$critical,coverage=inf$axis_coverage,
      old_exclusions=sum(old$contrasts$significant),new_exclusions=sum(inf$contrasts$significant),
      estimates_unchanged=TRUE,standard_errors_unchanged=TRUE,support_unchanged=TRUE)
    atomic_rds(inf,paste0('results/data/',stem,'_inference.rds'))
    allcontrasts[[axis]] <- mutate(inf$contrasts,axis=axis)
    if(label=='primary' && axis=='phase') {
      av <- interval_average_contrasts(inf,b)
      save_table(av,'primary_phase_interval_averages')
      wf <- bind_rows(lapply(seq_len(nrow(av)),function(i) {
        ix <- which(wrap24(inf$cell_design$x-av$start[i])<av$duration_hours[i]-1e-10)
        bind_rows(lapply(c('Work','Free'),function(context) {
          D <- inf$pattern_matrices[[paste0('Q4_',context)]]-inf$pattern_matrices[[paste0('Q1_',context)]]
          z <- linear_summary(matrix(colMeans(D[ix,,drop=FALSE]),nrow=1),coef(b$fit),covariance_info(b$fit)$V)
          bind_cols(tibble(interval=i,start=av$start[i],end=av$end[i],context=context),z) |>
            mutate(ratio=10^estimate,ratio_lower=10^lower,ratio_upper=10^upper)
        }))
      }))
      save_table(wf,'primary_phase_work_free_interval_averages')
      print(av |> select(start,end,ratio,ratio_lower,ratio_upper),width=Inf)
      print(wf |> select(start,end,context,ratio,ratio_lower,ratio_upper),width=Inf)
    }
    rm(b,old,inf);gc()
  }
  z <- bind_rows(allcontrasts);placement <- if(label=='primary')'Near-eye'else'Chest'
  for(context in c('Overall','Work','Free'))for(family in c('Evening minus morning','Extreme minus moderate'))
    save_plot(difference_figure(z,family,context,placement),paste0(label,'_',gsub(' ','_',tolower(family)),'_',tolower(context)),12,11)
  for(family in c('Free minus work','Interactions'))
    save_plot(difference_figure(z,family,placement=placement),paste0(label,'_',gsub(' ','_',tolower(family))),12,14)
}
save_table(bind_rows(checks),'revision7_inference_changes','results/provenance')
# Keep the compact manuscript-source exports consistent with the full tables.
phase_summary <- read.csv('results/tables/primary_phase_significant_windows.csv') |>
  filter(contrast=='evening_morning_Extreme_Overall') |>
  mutate(location=if_else(end<0,'Before published DLMO guide','After MSFsc'),
    evening_over_morning=mean_ratio,morning_over_evening=1/mean_ratio,
    mean_percent_deviation=100*(mean_ratio-1))
save_table(phase_summary,'summary_phase_intervals')
log_step('Refreshing noncyclic sensitivity with the same coordinate-specific coverage')
noncyclic_sz_outputs()
gc()
log_step('Refreshing MEQ inference on the four existing fitted models')
curves <- contrasts <- windows <- list();i <- 0L
meq_paths <- list.files('results/models',pattern='^meq_.*_(clock|phase)\\.rds$',full.names=TRUE)
stopifnot(length(meq_paths)==4L)
for(path in meq_paths) {
  i <- i+1L;b <- readRDS(path);inf <- meq_inference(b)
  curves[[i]] <- inf$curves;contrasts[[i]] <- inf$contrasts;windows[[i]] <- inf$windows
  rm(b,inf);gc()
}
curves <- bind_rows(curves);contrasts <- bind_rows(contrasts);windows <- bind_rows(windows)
save_table(curves,'meq_all_curves');save_table(contrasts,'meq_all_contrasts');save_table(windows,'meq_candidate_windows')
save_plot(meq_factor_figure(contrasts),'meq_factor_comparison',14,9)
log_step('Refreshing all photoperiod interval summaries, support counts and figures')
run_daylight(readRDS('results/data/prepared.rds'),validation=FALSE)
save_table(read.csv('results/tables/daylight_primary_interval_summaries.csv'),
  'summary_daylength_intervals')
after <- tools::md5sum(protected)
stopifnot(identical(before,after))
save_table(tibble(file=names(before),before_md5=unname(before),after_md5=unname(after),unchanged=before==after),
  'revision7_fitted_input_preservation','results/provenance')
writeLines(c('Command: Rscript --vanilla results/provenance/revision7_inference.R',
  paste('Executed:',Sys.time()),
  'Inputs: existing primary, chest, noncyclic-sz, MEQ, and photoperiod fitted models; prepared.rds; validation registry.',
  'Same fitted coefficients, covariance matrices, prediction distributions, observations, seeds and support rules.',
  '95% simultaneous coverage for the full contrast family within each coordinate; classification sensitivities evaluated separately.',
  'Photoperiod pointwise intervals remain 95%; interval summaries use the updated primary windows.',
  capture.output(sessionInfo())),'results/provenance/revision7_inference_runtime.txt')
log_step('Completed coverage revision; all fitted models and protected inputs unchanged')
