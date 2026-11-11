# Independent checks of the daylight model and analysis helpers.
# Run from the secondary-analysis project root:
# Rscript --vanilla results/provenance/daylight_independent_verify.R
# Updated 23 October 2026 for pointwise inference and the two primary solar windows.
# It reads cached models and writes only daylight_independent_* provenance files.
source('R/setup.R')
for (module in c('data','models','inference','validation','daylight')) {
  source(paste0('R/',module,'.R'))
}
command <- 'Rscript --vanilla results/provenance/daylight_independent_verify.R'
checks <- list()
check <- function(label, ok) {
  stopifnot(isTRUE(ok))
  checks[[length(checks)+1L]] <<- tibble(check=label,passed=TRUE)
}
prepared <- readRDS('results/data/prepared.rds')
raw <- prepared$frames$glasses
fold_registry <- read.csv('results/validation/participant_folds.csv')
fold_counts <- bind_rows(lapply(1:5,function(i) {
  assessment <- raw$participant_key %in%
    fold_registry$participant_key[fold_registry$fold==i]
  training <- raw[!assessment,]
  test <- raw[assessment,]
  prep <- fit_preprocess(training)
  original <- readRDS(paste0('results/validation/new_participant_fold_',i,'_phase.rds'))
  check(paste('Fold',i,'has disjoint participants'),
        length(intersect(training$participant_key,test$participant_key))==0L)
  check(paste('Fold',i,'training-only quartiles match primary validation'),
        isTRUE(all.equal(prep$cutoffs,original$cutoffs)))
  check(paste('Fold',i,'training-only photoperiod centers match primary validation'),
        isTRUE(all.equal(prep$centers,original$centers)))
  check(paste('Fold',i,'assessment sites occur in training'),
        all(unique(test$site)%in%unique(training$site)))
  tibble(fold=i,training_participants=n_distinct(training$participant_key),
         assessment_participants=n_distinct(test$participant_key))
}))

bundle <- readRDS('results/models/daylight_phase.rds')
d <- bundle$data
reference <- reference_days(d)
check('Reference weights sum to one',abs(sum(reference$w)-1)<1e-12)
check('Each participant has equal total reference weight',
      max(abs(tapply(reference$w,reference$participant_key,sum)-
              1/n_distinct(d$participant_key)))<1e-12)
probs <- c(.25,.5,.75)
targets <- daylight_weighted_quantile(reference$photo,reference$w,probs)
check('Scenario values are inverse weighted empirical quantiles',
      all(vapply(seq_along(targets),function(i) {
        sum(reference$w[reference$photo<targets[i]])<probs[i]+1e-12 &&
          sum(reference$w[reference$photo<=targets[i]])>=probs[i]-1e-12
      },logical(1))))
design <- list(targets=tibble(scenario=c('Shorter','Typical','Longer'),photo=targets),
               free_weight=sum(reference$w[reference$daytype=='Free']))
matrices <- daylight_matrices(bundle,design)
smooth <- Filter(function(s)identical(s$by,'photo'),bundle$fit$smooth)[[1]]
columns <- smooth$first.para:smooth$last.para
nd <- prediction_template(bundle,9)
nd$t <- seq(0,24,length.out=9)
nd$photo <- 0
zero <- population_matrix(bundle,nd)
nd$photo <- 1
one <- population_matrix(bundle,nd)
nd$photo <- 2
two <- population_matrix(bundle,nd)
check('Added cyclic sz term has 33 constrained columns',length(columns)==33L)
check('Added term is zero at centered photoperiod zero',
      max(abs(zero[,columns]))<1e-12)
check('Added term is linear in the numeric photoperiod modifier',
      max(abs(two[,columns]-2*one[,columns]))<1e-12)
quartile_matrices <- lapply(paste0('Q',1:4),function(q) {
  nd$quartile <- q
  population_matrix(bundle,nd)[,columns,drop=FALSE]
})
check('Added-term lpmatrix sums to zero over quartiles',
      max(abs(Reduce(`+`,quartile_matrices)))<1e-10)
check('Full population lpmatrix is cyclic at zero and 24 hours',
      max(abs(one[1,]-one[9,]))<1e-8)
other_columns <- setdiff(seq_len(ncol(matrices$D$Longer_minus_shorter)),columns)
check('Longer-minus-shorter contrast cancels every other model term',
      max(abs(matrices$D$Longer_minus_shorter[,other_columns]))<1e-10)

direct <- 0
for (day in c('Work','Free')) {
  nd <- prediction_template(bundle,length(matrices$x))
  nd$t <- wrap24(matrices$x)
  nd$photo <- design$targets$photo[1]
  nd$quartile <- 'Q4'
  nd$daytype <- day
  nd$direction <- 'Evening'
  nd$extremity <- 'Extreme'
  pp <- predict(bundle$fit,newdata=set_factors(nd,bundle),
                exclude=exclude_labels(bundle$fit),discrete=FALSE)
  weight <- if(day=='Free')design$free_weight else 1-design$free_weight
  direct <- direct+weight*pp
}
check('Direct predictions equal the weighted Shorter/Q4 design prediction',
      max(abs(direct-drop(matrices$X$Shorter_Q4%*%coef(bundle$fit))))<1e-8)

# Test window denominators without writing or replacing the author's tables.
window_function <- daylight_windows
environment(window_function) <- list2env(
  list(save_table=function(x,...)invisible(x)),
  parent=environment(daylight_windows))
x <- seq(-12,12-cfg$step,by=cfg$step)
for (wrap in c(FALSE,TRUE)) {
  synthetic <- tibble(contrast='Shorter',x=x,
    significant=if(wrap)x>=11.5|x<=(-11.5)else x>=1.333&x<=1.667,
    estimate=.1,ratio=10^.1,minimum_cell_participants=10L)
  synthetic$pointwise_significant<-synthetic$significant
  window <- window_function(list(data=d),list(contrasts=synthetic),design)
  label <- if(wrap)'Circular-boundary window'else'Ordinary window'
  check(paste(label,'has the intended boundary status'),
        nrow(window)==1L && window$crosses_boundary==wrap)
  for (q in c('Q1','Q4')) {
    rows <- d[as.character(d$quartile)==q &
      abs(d$photo-design$targets$photo[1])<=1 &
      wrap24(d$phase_display-window$start_phase)<window$duration_hours-1e-10,]
    nr <- count_sample(rows)
    check(paste(label,q,'counts match exact observed phase coordinates'),
      nr$participants==window[[paste0('participants_Shorter_',q)]] &&
      nr$participant_days==window[[paste0('participant_days_Shorter_',q)]] &&
      nr$half_hours==window[[paste0('half_hours_Shorter_',q)]])
  }
}
site <- d |> filter(quartile%in%c('Q1','Q4')) |>
  summarise(participants=n_distinct(participant_key),.by=c(site,quartile)) |>
  complete(site,quartile=c('Q1','Q4'),fill=list(participants=0L))
check('No observed site has at least five participants in both extreme groups',
      !any((site |> summarise(both=all(participants>=5),.by=site))$both))

# Independently evaluate the new reported bounds using saved model covariance.
reported<-read.csv('results/tables/daylight_contrasts.csv')
V<-covariance_info(bundle$fit)$V;beta<-coef(bundle$fit)
for(id in names(matrices$D)) {
  X<-matrices$D[[id]];eta<-drop(X%*%beta)
  se<-sqrt(pmax(0,rowSums((X%*%V)*X)))
  z<-reported[reported$contrast==id,];z<-z[order(z$x),]
  check(paste(id,'pointwise bounds equal direct coefficient-covariance calculations'),
    max(abs(z$pointwise_lower-(eta-qnorm(.975)*se)))<1e-9 &&
      max(abs(z$pointwise_upper-(eta+qnorm(.975)*se)))<1e-9)
  check(paste(id,'pointwise exclusions equal supported zero exclusions'),
    identical(z$pointwise_significant,z$supported & (z$pointwise_lower>0|z$pointwise_upper<0)))
}
moderation<-read.csv('results/tables/daylight_pointwise_moderation_exclusions.csv')
check('Moderation export contains the same complete 144-point grid',
      nrow(moderation)==144L && max(abs(moderation$x-matrices$x))<1e-12)
averages<-read.csv('results/tables/daylight_primary_interval_summaries.csv')
for(i in seq_len(nrow(averages))) {
  z<-averages[i,];ix<-which(wrap24(matrices$x-z$start_phase)<z$duration_hours-1e-10)
  X<-colMeans(matrices$D[[z$contrast]][ix,,drop=FALSE])
  eta<-sum(X*beta);se<-sqrt(drop(crossprod(X,V%*%X)))
  check(paste('Interval',z$interval,z$contrast,'uses the covariance of the averaged design'),
    abs(z$ratio-10^eta)<1e-8 && abs(z$ratio_lower-10^(eta-qnorm(.975)*se))<1e-8 &&
      abs(z$ratio_upper-10^(eta+qnorm(.975)*se))<1e-8)
}

# Rebuild solar availability directly as a participant-day by phase matrix.
context<-readRDS(input_paths()[['glasses_context']])
solar_days<-context |> distinct(site,Id,position,local_date,civil_dawn_wall_minute,civil_dusk_wall_minute)
keys<-c('site','Id','position','local_date')
check('Solar metadata has one record per participant local date',!anyDuplicated(solar_days[keys]))
days<-d |> distinct(participant_key,participant_day_key,site,Id,position,local_date,quartile,msf_hours,photoperiod_hours) |>
  left_join(solar_days,by=keys)
check('Every retained day has complete solar metadata',
  nrow(days)==n_distinct(d$participant_day_key)&&!anyNA(days$civil_dawn_wall_minute)&&!anyNA(days$civil_dusk_wall_minute))
weights<-1/as.numeric(table(days$participant_key)[days$participant_key]);o<-order(days$photoperiod_hours)
median_day<-days$photoperiod_hours[o[which(cumsum(weights[o])/sum(weights)>=.5)[1]]]
check('Absolute-daylength split is the participant-weighted empirical median',
  abs(median_day-read.csv('results/tables/daylight_civil_definition.csv')$weighted_median_hours)<1e-12)
days$daylength<-ifelse(days$photoperiod_hours<=median_day,'Shorter','Longer')
summaries<-read.csv('results/tables/daylight_civil_primary_intervals.csv')
primary_windows<-read.csv('results/tables/primary_phase_significant_windows.csv') |>
  filter(contrast=='evening_morning_Extreme_Overall') |> arrange(start)
check('Civil summaries have exactly the two primary windows and four cells per window',
  nrow(summaries)==8L && all(table(summaries$interval)==4L) &&
    max(abs(summaries$start_phase-primary_windows$start[summaries$interval]))<1e-12 &&
    max(abs(summaries$end_phase-primary_windows$end[summaries$interval]))<1e-12)
for(i in seq_len(nrow(summaries))) {
  z<-summaries[i,];dd<-days[as.character(days$quartile)==z$quartile & days$daylength==z$daylength,]
  xx<-matrices$x[wrap24(matrices$x-z$start_phase)<z$duration_hours-1e-10]
  clock<-outer(dd$msf_hours,xx,function(a,b)((a+b)%%24)*60)
  since_dawn<-sweep(clock,1,dd$civil_dawn_wall_minute,'-')%%1440
  duration<-(dd$civil_dusk_wall_minute-dd$civil_dawn_wall_minute)%%1440
  available<-sweep(since_dawn,1,duration,'<=')
  ww<-1/as.numeric(table(dd$participant_key)[dd$participant_key]);ww<-ww/sum(ww)
  estimate<-sum(ww*rowMeans(available))
  check(paste('Solar interval',z$interval,z$quartile,z$daylength,'matches direct weighted day-phase calculations'),
    abs(estimate-z$mean_availability)<1e-12 &&
      z$participant_days==nrow(dd) && z$participants==n_distinct(dd$participant_key) && z$grid_points==length(xx))
  observed<-d[d$participant_day_key %in% dd$participant_day_key & wrap24(d$phase_display-z$start_phase)<z$duration_hours-1e-10,]
  check(paste('Solar interval',z$interval,z$quartile,z$daylength,'separately records observed support'),
    z$interval_observed_participants==n_distinct(observed$participant_key)&&
      z$interval_observed_days==n_distinct(observed$participant_day_key)&&z$interval_observed_half_hours==nrow(observed))
}

input_files <- c('results/data/prepared.rds','results/models/daylight_phase.rds',
  'results/validation/participant_folds.csv',
  paste0('results/validation/new_participant_fold_',1:5,'_phase.rds'),
  paste0('R/',c('setup','data','models','inference','validation','daylight'),'.R'),
  'results/tables/daylight_contrasts.csv','results/tables/daylight_primary_interval_summaries.csv',
  'results/tables/daylight_civil_primary_intervals.csv','results/tables/daylight_civil_definition.csv',
  'results/tables/primary_phase_significant_windows.csv',input_paths()[['glasses_context']],
  'results/provenance/daylight_independent_verify.R')
write.csv(tibble(file=input_files,md5=unname(tools::md5sum(input_files))),
          'results/provenance/daylight_independent_inputs.csv',row.names=FALSE)
write.csv(bind_rows(checks),'results/provenance/daylight_independent_checks.csv',row.names=FALSE)
write.csv(fold_counts,'results/provenance/daylight_independent_folds.csv',row.names=FALSE)
writeLines(c(paste('Command:',command),paste('Working directory:',getwd()),
  paste('Executed:',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')),
  paste('OpenMP preload:',omp),
  paste('Passing checks:',length(checks)),
  'Scope: cached full-model verification only; no model fitting or coefficient draws.',
  '',capture.output(sessionInfo())),
  'results/provenance/daylight_independent_runtime.txt')
cat(length(checks),'independent checks passed.\n')
