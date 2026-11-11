# R-only scientific checks of revision 5. Run from the project root.
source('R/setup.R')
for(m in c('models','inference','plots','daylight','phase_response')) source(paste0('R/',m,'.R'))
checks<-list()
check<-function(name,passed,error=NA_real_) {
  checks[[length(checks)+1L]]<<-tibble(check=name,passed=isTRUE(passed),error=error)
  if(!isTRUE(passed))stop(name)
}
near<-function(name,a,b,tol=1e-10) {error<-max(abs(a-b));check(name,all(is.finite(a))&&all(is.finite(b))&&error<tol,error)}
# The source coefficients, published-coordinate import and original zero-light
# results provide regression references independent of the new centering code.
old<-readRDS('results/provenance/revision5-original/results/data/phase_response_results.rds')
z<-readRDS('results/data/phase_response_results.rds')
b<-readRDS('results/models/primary_phase.rds')
i<-readRDS('results/data/primary_phase_inference.rds')
near('Published PRC coefficients unchanged',z$prc$coefficients,old$prc$coefficients)
keys<-c('quartile','context','dlmo_offset','prc_baseline','exposure_scale')
a<-z$summary |> filter(exposure_scale=='Zero-exposure centered') |> arrange(across(all_of(keys)))
c<-old$summary |> filter(quartile %in% c('Q1','Q4'),exposure_scale=='Zero-exposure centered') |> arrange(across(all_of(keys)))
check('Zero-exposure regression keys identical',identical(a[keys],c[keys]))
near('Zero-exposure estimates and intervals reproduce the previous analysis',
  as.matrix(a[c('estimate','se','lower','upper')]),as.matrix(c[c('estimate','se','lower','upper')]))
check('Only extreme groups with zero-exposure centering remain',
  setequal(z$summary$quartile,c('Q1','Q4')) && identical(unique(z$summary$exposure_scale),'Zero-exposure centered') &&
  nrow(z$summary)==60L && nrow(z$contrasts)==30L)
# Independently evaluate all integrals and variances on the full saved design.
beta<-coef(b$fit);V<-if(!is.null(b$fit$Vc))b$fit$Vc else b$fit$Vp
x<-i$cell_design$x;A<-list();pred<-lapply(i$pattern_matrices,function(X)drop(X%*%beta))
for(n in seq_len(nrow(z$summary))) {
  row<-z$summary[n,];key<-row$key
  kernel<-drop(z$prc$basis(x-row$dlmo_offset)%*%z$prc$coefficients)+if(row$prc_baseline=='Drift corrected').54 else 0
  weights<-kernel/sum(abs(kernel));baseline<- -1
  estimate<-sum(weights*(pred[[key]]-baseline))
  aa<-colSums(i$pattern_matrices[[key]]*weights)
  se<-sqrt(drop(aa%*%V%*%aa));ci<-estimate+c(-1,1)*qnorm(.975)*se
  near(paste('Full-design integral and covariance',n),c(estimate,se,ci),unlist(row[c('estimate','se','lower','upper')]));
  A[[paste(row$dlmo_offset,row$prc_baseline,key,sep='|')]]<-aa
}
for(n in seq_len(nrow(z$contrasts))) {
  row<-z$contrasts[n,];kk<-paste(c('Q1','Q4'),row$context,sep='_')
  selected<-z$summary |> filter(dlmo_offset==row$dlmo_offset,prc_baseline==row$prc_baseline,
    context==row$context,exposure_scale==row$exposure_scale) |> arrange(quartile)
  estimate<-diff(selected$estimate)
  aa<-A[[paste(row$dlmo_offset,row$prc_baseline,kk[2],sep='|')]]-A[[paste(row$dlmo_offset,row$prc_baseline,kk[1],sep='|')]]
  se<-sqrt(drop(aa%*%V%*%aa));ci<-estimate+c(-1,1)*qnorm(.975)*se
  near(paste('Joint group-contrast covariance',n),c(estimate,se,ci),unlist(row[c('estimate','se','lower','upper')]))
}
check('Drift-corrected central indices are advance-weighted',
  all((z$summary |> filter(dlmo_offset== -6.45,prc_baseline=='Drift corrected'))$lower>0))
raw<-z$summary |> filter(dlmo_offset== -6.45,prc_baseline=='Published raw')
check('Every published-raw central index interval includes zero',all(raw$lower<0 & raw$upper>0))
check('Obsolete median-centered results and fields are absent',
  !any(c('daily_medians','median_references') %in% names(z)) &&
  !file.exists('results/tables/phase_response_median_references.csv') &&
  !file.exists('results/data/phase_response_daily_medians.rds'))
check('Every central-alignment evening-minus-morning point estimate is positive',
  all((z$contrasts |> filter(dlmo_offset== -6.45))$estimate>0))
check('Every central-alignment group contrast includes zero',
  all((z$contrasts |> filter(dlmo_offset== -6.45))$lower<0) && all((z$contrasts |> filter(dlmo_offset== -6.45))$upper>0))
exception<-z$contrasts |> filter(reference_alignment,lower>0|upper<0)
check('Exactly the reported early-guide raw-PRC work-day exception',nrow(exception)==1L &&
  exception$dlmo_offset== -6.95 && exception$context=='Work' &&
  exception$prc_baseline=='Published raw' && exception$exposure_scale=='Zero-exposure centered')
save_table(exception,'revision5_prc_alignment_exception','results/provenance')
rm(b,i,old,V,A);gc()
# Daylength estimates, targets and support must remain identical to the reviewed
# values; this revision adds explanatory annotations only.
for(spec in c('fixed','quartiles','p10_p90')) {
  prefix<-if(spec=='fixed')'daylight'else paste0('daylight_',spec)
  des<-readRDS(paste0('results/data/',prefix,'_design.rds'))
  inf<-readRDS(paste0('results/data/',prefix,'_inference.rds'))
  check(paste('Every phase point meets local n=4:',spec),all(inf$contrasts$supported) && min(inf$contrasts$minimum_cell_participants)>=4L)
  check(paste('Pointwise significance matches intervals:',spec),
    identical(inf$contrasts$pointwise_significant,inf$contrasts$supported &
      (inf$contrasts$pointwise_lower>0|inf$contrasts$pointwise_upper<0)))
  mod<-inf$contrasts |> filter(contrast=='Longer_minus_shorter')
  check(paste('Thirteen significant phase-grid coordinates:',spec),sum(mod$pointwise_significant)==13L)
  near(paste('Reported significant phase range:',spec),range(mod$x[mod$pointwise_significant]),c(7/6,19/6))
  if(spec=='fixed')near('Main targets are exactly minus one, zero and plus one',des$targets$photo,c(-1,0,1))
  if(spec=='p10_p90')check('Four-person cells included in P10/P90',min(inf$contrasts$minimum_cell_participants)==4L)
  if(spec=='quartiles') {
    old<-readRDS('results/provenance/revision4-original/results/data/daylight_inference.rds')
    a<-inf$contrasts |> arrange(contrast,x);c<-old$contrasts |> arrange(contrast,x)
    near('Original quartile contrast estimates, SEs and pointwise bounds preserved',
      as.matrix(a[c('estimate','se','pointwise_lower','pointwise_upper')]),as.matrix(c[c('estimate','se','pointwise_lower','pointwise_upper')]))
  }
}
check('Primary support threshold remains five',cfg$support_n==5L)
files<-list.files('results/tables',pattern='^daylight.*[.]csv$',full.names=TRUE)
for(f in files) {
  old_path<-file.path('results/provenance/revision5-original',f)
  check(paste('Daylength numerical table unchanged:',basename(f)),
    file.exists(old_path) && identical(read.csv(f),read.csv(old_path)))
}
# Check phase locations against the primary contrasts, independent of plot styling.
primary<-read.csv('results/tables/primary_phase_significant_windows.csv')
primary<-primary[primary$contrast=='evening_morning_Extreme_Overall',]
primary<-primary[order(primary$start),]
band<-primary_phase_interval_band()$geom_params
near('Primary interval overlay starts',band$starts,primary$start)
near('Primary interval overlay ends',band$ends,primary$end)
near('DLMO guide remains separate',c(dlmo_reference_band()$geom_params$xmin,
  dlmo_reference_band()$geom_params$xmax),c(-6.95,-5.95))
moderation<-read.csv('results/tables/daylight_pointwise_contrasts.csv')
pre<-moderation[moderation$contrast=='Longer_minus_shorter' & moderation$x< -6.95,]
trough<-pre[which.min(pre$ratio),]
near('Pre-DLMO trough timing',trough$x,-25/3)
check('Pre-DLMO trough is later than the primary early interval',trough$x>primary$end[1])
check('Pre-DLMO trough interval still includes one',trough$pointwise_ratio_lower<1 && trough$pointwise_ratio_upper>1)
check('Pre-DLMO reported rounded ratio and bounds',
  identical(sprintf('%.2f',unlist(trough[c('ratio','pointwise_ratio_lower','pointwise_ratio_upper')])),c('0.41','0.16','1.04')))
save_table(trough,'revision5_predlmo_trough','results/provenance')

# Verify equal-participant photoperiod percentiles using an independent weighted
# empirical CDF, and proportional rescaling of the moderation estimates/SEs.
b<-readRDS('results/models/daylight_phase.rds');days<-unique(b$data[c('participant_key','participant_day_key','photo')])
n<-table(days$participant_key);w<-1/as.numeric(n[days$participant_key])/length(n);o<-order(days$photo);cdf<-cumsum(w[o])/sum(w)
quant<-function(p)days$photo[o[which(cdf>=p)[1]]]
targets<-read.csv('results/tables/daylight_sensitivity_targets.csv')
for(spec in c('quartiles','p10_p90')) {
  a<-targets[targets$specification==spec,]
  near(paste('Independent weighted percentile targets:',spec),a$photo,vapply(a$probability,quant,numeric(1)))
}
base<-read.csv('results/tables/daylight_pointwise_contrasts.csv');base<-base[base$contrast=='Longer_minus_shorter',]
for(spec in c('quartiles','p10_p90')) {
  t<-targets[targets$specification==spec,];span<-diff(range(t$photo))/2
  a<-read.csv(paste0('results/tables/daylight_',spec,'_pointwise_contrasts.csv'));a<-a[a$contrast=='Longer_minus_shorter',]
  near(paste('Linear moderation rescaling:',spec),c(a$estimate,a$se),c(base$estimate,base$se)*span)
}
rm(b);gc()
# All fitted model and primary inference, subset and validation files must be
# byte-identical to the pre-revision manifest.
before<-read.csv('results/provenance/revision5_models_before.csv')
after<-transform(before,md5=unname(tools::md5sum(path)))
check('All protected fitted models and primary analysis caches unchanged',identical(before,after))
save_table(after,'revision5_models_after','results/provenance')
save_table(bind_rows(checks),'revision5_analysis_checks','results/provenance')
writeLines(c('Command: Rscript --vanilla results/provenance/revision5_verify_analysis.R',
  paste('Completed:',Sys.time()),paste('Passed:',length(checks)),capture.output(sessionInfo())),
  'results/provenance/revision5_analysis_runtime.txt')
cat(length(checks),'analytical checks passed.\n')
