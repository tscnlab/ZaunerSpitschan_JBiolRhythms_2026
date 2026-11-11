# Read-only scientific audit. Only revision_floor_audit* provenance files are written.
# Command: Rscript --vanilla results/provenance/revision_floor_audit.R
source('R/setup.R')
source('R/models.R')
source('R/inference.R')
source('R/daylight.R')
prefix <- 'results/provenance/revision_floor_audit'
write_audit <- function(x, suffix) write.csv(x, paste0(prefix, '_', suffix, '.csv'), row.names=FALSE)
inputs <- c('R/setup.R','R/models.R','R/inference.R','R/daylight.R','R/sensitivity.R','R/attribution.R',
  paste0('results/models/',c('primary_clock','primary_phase','variance_weighted_clock','variance_weighted_phase','daylight_phase'),'.rds'),
  paste0('results/subsets/',c('clock','phase'),'_allocation.rds'),
  'results/data/daylight_design.rds','results/tables/model_fit.csv',
  'results/tables/daylight_primary_interval_summaries.csv',
  'results/tables/primary_phase_significant_windows.csv','results/tables/published_dlmo_reference.csv',
  'results/tables/primary_clock_individual_timing.csv','results/tables/primary_phase_individual_timing.csv',
  'results/tables/sensitivity_summary.csv')
write_audit(data.frame(path=inputs,sha256=vapply(inputs,function(p)digest::digest(file=p,algo='sha256'),character(1))), 'inputs')
writeLines(c('Command: Rscript --vanilla results/provenance/revision_floor_audit.R',
  'No models were refitted. Existing saved fits and inference are read-only.',
  paste('OpenMP preload:',omp),capture.output(sessionInfo()),
  paste('mgcv',packageVersion('mgcv'),'dplyr',packageVersion('dplyr'),'gratia',packageVersion('gratia'))),paste0(prefix,'_session.txt'))
fits <- spreads <- summaries <- counts <- fitstats <- record_hashes <- list()
for(axis in c('clock','phase')) for(kind in c('primary','variance_weighted')) {
  obj <- readRDS(paste0('results/models/',kind,'_',axis,'.rds'))
  b <- if(kind=='primary')obj else obj$bundle
  f<-b$fit;d<-b$data;n<-nrow(d);mu<-as.numeric(fitted(f));e<-as.numeric(residuals(f,type='response'));z<-as.numeric(f$std.rsd)
  stopifnot(length(z)==n,max(abs(d$response-log10(d$metric_value_lx+.1)))<1e-12,
    max(abs(e-(d$response-mu)))<1e-10,all(d$metric_value_lx>=0))
  previous<-c(NA_real_,head(d$metric_value_lx,-1))
  d<-mutate(d,fit_kind=kind,axis=axis,mu=mu,residual=e,innovation=z,
    previous_exposure=previous,linked=!AR_start,
    fitted_decile=ntile(mu,10),observed_zero=metric_value_lx==0,
    observed_below_01=metric_value_lx<.1,observed_below_1=metric_value_lx<1)
  id<-paste(kind,axis,sep='_')
  record_hashes[[id]]<-digest::digest(d[,c('participant_key','participant_day_key','clock_bin','source_utc_start','response')])
  counts[[id]]<-d |> summarise(fit_kind=.env$kind,axis=.env$axis,
    participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),
    zero_half_hours=sum(observed_zero),zero_share=mean(observed_zero),
    below_01_half_hours=sum(observed_below_01),below_01_share=mean(observed_below_01),
    below_1_half_hours=sum(observed_below_1),below_1_share=mean(observed_below_1),
    fitted_below_floor=sum(mu< -1),minimum_fitted=min(mu),
    floor_relation_max_error=max(abs(residual[observed_zero]-(-1-mu[observed_zero]))))
  fits[[id]]<-d |> summarise(mean_fitted=mean(mu),minimum_fitted=min(mu),maximum_fitted=max(mu),
    zero_share=mean(observed_zero),below_01_share=mean(observed_below_01),below_1_share=mean(observed_below_1),
    raw_sd=sd(residual),raw_IQR=IQR(residual),innovation_sd=sd(innovation),
    participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),
    .by=c(fit_kind,axis,fitted_decile))
  definitions <- list(all=rep(TRUE,n),positive=d$metric_value_lx>0,
    at_least_01_lx=d$metric_value_lx>=.1,at_least_1_lx=d$metric_value_lx>=1,
    at_least_10_lx=d$metric_value_lx>=10,
    positive_current_previous=d$linked & d$metric_value_lx>0 & previous>0,
    at_least_1_lx_current_previous=d$linked & d$metric_value_lx>=1 & previous>=1,
    fitted_at_least_1_lx=mu>=log10(1.1))
  for(label in names(definitions)) {
    ds<-d[which(definitions[[label]]),,drop=FALSE]
    # Re-bin within each retained stratum so its exposure range has equal support.
    ss<-ds |> mutate(stratum=label,stratum_decile=ntile(mu,10)) |>
      summarise(mean_fitted=mean(mu),minimum_fitted=min(mu),maximum_fitted=max(mu),
        raw_sd=sd(residual),raw_IQR=IQR(residual),innovation_sd=sd(innovation),
        participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),
        .by=c(fit_kind,axis,stratum,stratum_decile))
    stopifnot(all(ss$half_hours>=100),all(ss$participants>=5))
    spreads[[paste(id,label)]]<-ss
    summaries[[paste(id,label)]]<-data.frame(fit_kind=kind,axis=axis,stratum=label,
      participants=n_distinct(ds$participant_key),participant_days=n_distinct(ds$participant_day_key),half_hours=nrow(ds),
      raw_sd_min=min(ss$raw_sd),raw_sd_max=max(ss$raw_sd),raw_sd_ratio=max(ss$raw_sd)/min(ss$raw_sd),
      raw_IQR_ratio=max(ss$raw_IQR)/min(ss$raw_IQR),innovation_sd_min=min(ss$innovation_sd),
      innovation_sd_max=max(ss$innovation_sd),innovation_sd_ratio=max(ss$innovation_sd)/min(ss$innovation_sd))
  }
  if(kind=='primary') {
    ordinary<-1-sum(e^2)/sum((d$response-mean(d$response))^2)
    # mgcv summary.gam definition, calculated directly to avoid repeated smooth tests.
    adjusted<-1-var(e)*((n-1)/(n-sum(f$edf)))/var(d$response)
    fitstats[[axis]]<-data.frame(axis=axis,ordinary_R2=ordinary,adjusted_R2=adjusted,
      deviance_explained=1-f$deviance/f$null.deviance,
      participants=n_distinct(d$participant_key),participant_days=n_distinct(d$participant_day_key),half_hours=n)
    tab<-read.csv('results/tables/model_fit.csv') |> filter(model=='primary',.data$axis==.env$axis)
    stopifnot(abs(tab$ordinary_R2-ordinary)<1e-10,abs(tab$adjusted_R2-adjusted)<1e-10)
  }
  cat('Completed',id,'\n');flush.console();rm(obj,b,f,d);gc()
}
stopifnot(length(unique(unlist(record_hashes)))==1L)
write_audit(data.frame(model=names(record_hashes),identical_record_hash=unlist(record_hashes)), 'sample_agreement')
write_audit(bind_rows(counts),'exposure_shares')
write_audit(bind_rows(fits),'original_deciles')
write_audit(bind_rows(spreads),'stratum_deciles')
write_audit(bind_rows(summaries),'spread_summaries')
write_audit(bind_rows(fitstats),'R2')

# Joint bootstrap sums retain covariance. Never add marginal CI endpoints.
allocations<-list()
for(axis in c('clock','phase')) {
  a<-readRDS(paste0('results/subsets/',axis,'_allocation.rds'))
  stopifnot(max(abs(rowSums(a$bootstrap)-(a$full_bootstrap-a$baseline_bootstrap)))<1e-8,
    abs(sum(a$table$allocation)-(a$table$full_R2[1]-a$table$baseline_R2[1]))<1e-9)
  groupings<-list(chronotype=c('direction','extremity'),
    individual_variation=c('participant','participant_day'),all_seven=a$table$component)
  for(label in names(groupings)) {
    cols<-match(groupings[[label]],a$table$component);stopifnot(!anyNA(cols))
    draws<-rowSums(a$bootstrap[,cols,drop=FALSE]);ci<-quantile(draws,c(.025,.975),type=8)
    allocations[[paste(axis,label)]]<-data.frame(axis=axis,combined_group=label,
      allocation=sum(a$table$allocation[cols]),lower=ci[1],upper=ci[2],replicates=length(draws),
      interval='95% percentile interval from the joint site-participant-day fixed-prediction bootstrap; conditional on fitted models',
      participants=a$table$participants[1],participant_days=a$table$participant_days[1],half_hours=a$table$half_hours[1])
  }
}
write_audit(bind_rows(allocations),'combined_shapley')

timing<-bind_rows(lapply(c('clock','phase'),function(axis) {
  path<-paste0('results/tables/primary_',axis,'_individual_timing.csv');z<-read.csv(path)
  data.frame(axis=axis,participants=nrow(z),localized_peaks=sum(z$localized),
    unlocalized_peaks=sum(!z$localized),participant_days=sum(z$participant_days),half_hours=sum(z$half_hours),
    median_peak_interval_width=median(z$peak_interval_width),
    minimum_peak_interval_width=min(z$peak_interval_width),maximum_peak_interval_width=max(z$peak_interval_width))
}))
write_audit(timing,'peak_localization')
deletion<-read.csv('results/tables/sensitivity_summary.csv') |> filter(sensitivity=='delete_site')
write_audit(deletion,'site_deletion')

guide<-read.csv('results/tables/published_dlmo_reference.csv')
win<-read.csv('results/tables/primary_phase_significant_windows.csv') |>
  filter(contrast=='evening_morning_Extreme_Overall',end<0)
stopifnot(nrow(win)==1)
distance<-data.frame(window_start=win$start,window_end=win$end,
  guide_start=min(guide$relative_DLMO),guide_end=max(guide$relative_DLMO),
  closest_gap_hours=min(guide$relative_DLMO)-win$end,
  furthest_gap_hours=max(guide$relative_DLMO)-win$start,
  window_duration_hours=win$duration_hours)
write_audit(distance,'early_DLMO_distance')

# Current conservative envelopes and an optional conditional interval-average CI.
# The selected windows come from the primary fit, so the latter does not account
# for interval selection and must not replace full-family simultaneous inference.
old<-read.csv('results/tables/daylight_primary_interval_summaries.csv')
if('simultaneous_envelope_lower' %in% names(old)) {
  old$ratio_lower<-old$simultaneous_envelope_lower
  old$ratio_upper<-old$simultaneous_envelope_upper
}
old<-old |> mutate(percent_deviation=100*(ratio-1),percent_lower_envelope=100*(ratio_lower-1),
  percent_upper_envelope=100*(ratio_upper-1),interval_type='Exponentiated mean of simultaneous log-contrast envelopes')
write_audit(old,'daylight_percent_envelopes')
b<-readRDS('results/models/daylight_phase.rds');design<-readRDS('results/data/daylight_design.rds')
mm<-daylight_matrices(b,design);vc<-covariance_info(b$fit);beta<-coef(b$fit)
av<-list()
for(i in seq_len(nrow(old))) {
  z<-old[i,];ix<-which(wrap24(mm$x-z$start_phase)<z$duration_hours-1e-10)
  X<-matrix(colMeans(mm$D[[z$contrast]][ix,,drop=FALSE]),nrow=1)
  ss<-linear_summary(X,beta,vc$V)
  stopifnot(abs(ss$estimate-z$mean_log10_difference)<1e-9)
  av[[i]]<-data.frame(interval=z$interval,start_phase=z$start_phase,end_phase=z$end_phase,
    contrast=z$contrast,grid_points=length(ix),estimate=ss$estimate,se=ss$se,
    ratio=10^ss$estimate,ratio_lower=10^ss$lower,ratio_upper=10^ss$upper,
    percent_deviation=100*(10^ss$estimate-1),percent_lower=100*(10^ss$lower-1),percent_upper=100*(10^ss$upper-1),
    covariance=vc$label,
    interval_type='95% conditional Gaussian interval for the average log contrast; selected interval treated as fixed; no multiplicity or selection adjustment')
}
write_audit(bind_rows(av),'daylight_interval_average')
print(bind_rows(counts));print(bind_rows(summaries));print(bind_rows(fitstats));print(bind_rows(allocations));print(distance)
print(bind_rows(av));cat('Audit complete. No fitted models or shared analysis files were changed.\n')
model_paths<-inputs[grepl('results/models/',inputs,fixed=TRUE)]
before<-read.csv(paste0(prefix,'_inputs.csv'))
unchanged<-data.frame(path=model_paths,sha256_after=vapply(model_paths,function(p)digest::digest(file=p,algo='sha256'),character(1)))
unchanged$unchanged<-unchanged$sha256_after==before$sha256[match(unchanged$path,before$path)]
stopifnot(all(unchanged$unchanged));write_audit(unchanged,'model_preservation')
