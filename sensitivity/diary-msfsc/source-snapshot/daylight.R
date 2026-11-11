# Exploratory photoperiod moderation. All artifacts have a daylight prefix;
# the primary models, estimands, data cache and validation registry are read-only.
daylight_version <- '1.1.0'
# Reporting changes are separate from the fitting signature. No refit is needed
# when the same covariance is summarized with pointwise or interval-average CIs.
daylight_reporting_version <- '3.0.0'
daylight_support_n <- 4L
# mgcv 1.9-4 otherwise discards the sz constraint for a varying numeric by.
# Keep its original constructor, basis, penalties and prediction method; force
# its sum-to-zero constraint to apply before multiplying by photoperiod.
smooth.construct.szphoto.smooth.spec <- function(object,data,knots) {
  out<-mgcv:::smooth.construct.sz.smooth.spec(object,data,knots)
  attr(out$C,'always.apply')<-TRUE
  out
}
# Retain the native sz declaration used by bam's discrete tensor bookkeeping.
smooth.info.szphoto.smooth.spec <- function(object) {
  object$tensor.possible<-TRUE
  object
}
daylight_formula <- function() update(model_formula(), . ~ . +
  s(t,quartile,by=photo,bs='szphoto',xt=list(bs='cc'),k=12))

daylight_weighted_quantile <- function(x,w,probs) {
  stopifnot(all(is.finite(x)),all(w>=0),sum(w)>0)
  o<-order(x); cdf<-cumsum(w[o])/sum(w)
  vapply(probs,function(p)x[o[which(cdf>=p)[1]]],numeric(1))
}

daylight_design <- function(d,specification=c('fixed','quartiles','p10_p90'),prefix='daylight') {
  specification<-match.arg(specification);ref<-reference_days(d)
  probs<-switch(specification,fixed=c(NA_real_,NA_real_,NA_real_),quartiles=c(.25,.5,.75),p10_p90=c(.1,.5,.9))
  targets<-tibble(scenario=c('Shorter','Typical','Longer'),probability=probs,
    photo=if(specification=='fixed')c(-1,0,1)else daylight_weighted_quantile(ref$photo,ref$w,probs))
  support<-bind_rows(lapply(seq_len(nrow(targets)),function(i) {
    d |> filter(abs(photo-targets$photo[i])<=1,quartile %in% c('Q1','Q4')) |>
      summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(quartile,daytype)) |>
      mutate(scenario=targets$scenario[i],photo_target=targets$photo[i],photo_tolerance_hours=1)
  }))
  stopifnot(nrow(support)==12L,all(support$participants>=daylight_support_n))
  site<-d |> filter(quartile %in% c('Q1','Q4')) |> mutate(quartile=as.character(quartile)) |>
    summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(site,quartile)) |>
    complete(site,quartile=c('Q1','Q4'),fill=list(participants=0L,participant_days=0L,half_hours=0L))
  save_table(targets,paste0(prefix,'_scenario_targets'));save_table(support,paste0(prefix,'_scenario_samples'))
  save_table(site,paste0(prefix,'_site_composition'))
  design<-list(targets=targets,free_weight=sum(ref$w[ref$daytype=='Free']),support=support,site=site,
    formula=deparse(daylight_formula()),version=daylight_version,reporting_version=daylight_reporting_version,
    specification=specification,prefix=prefix,support_n=daylight_support_n,
    inference_family=c(paste('Q4 minus Q1 at',targets$scenario),'Longer minus shorter of Q4 minus Q1'),
    local_support='At least four participants in each extreme-quartile by day-type cell, within one hour of the centered photoperiod target and the half-hour phase bin.',
    quantiles='Inverse weighted empirical CDF; equal participant weights, then equal weights for their observed days.',
    interpretation='Exploratory pooled association with site adjustment; targets and support threshold revised after examining the fit.')
  atomic_rds(design,paste0('results/data/',prefix,'_design.rds'))
  writeLines(c(paste('Exploratory reporting specification:',specification),
    paste('Formula:',paste(design$formula,collapse=' ')),
    paste('Targets:',paste(targets$scenario,format(targets$photo,digits=8),collapse='; ')),
    paste('Pooled free-day weight:',format(design$free_weight,digits=9)),design$quantiles,
    paste('Inference family:',paste(design$inference_family,collapse='; ')),
    'Pointwise 95% Gaussian bands and full-covariance interval-average CIs; exploratory and unadjusted.',
    design$local_support,design$interpretation,
    'Moderation is linear in centered photoperiod and shared across work and free days. Primary fits are unchanged.',
    'szphoto retains the original mgcv sz constructor, prediction method and sum-to-zero constraint.'),
    paste0('results/provenance/',prefix,'_design.txt'))
  design
}

daylight_fit <- function(raw,prep,label='daylight',rho_override=NULL,k=cfg$k,kp=cfg$kp) {
  d<-apply_preprocess(raw,prep,'phase')
  formula<-update(model_formula(k=k,kp=kp),. ~ . + s(t,quartile,by=photo,bs='szphoto',xt=list(bs='cc'),k=12))
  signature<-digest::digest(list(d,formula=deparse(formula),rho_override,daylight_version,mgcv=as.character(packageVersion('mgcv'))))
  path<-paste0('results/models/',label,'_phase.rds')
  if(file.exists(path)) {old<-readRDS(path);if(identical(old$signature,signature))return(old);rm(old);gc()}
  preliminary<-fit_bam(d,formula,0,paste(label,'phase preliminary'))
  set.seed(cfg$seed+700L)
  kc<-mgcv::k.check(preliminary,subsample=5000,n.rep=100)
  patterns<-d |> mutate(residual=residuals(preliminary),hour=floor(t)) |>
    summarise(mean_residual=mean(residual),sd_residual=sd(residual),participants=n_distinct(participant_key),
      participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(hour,quartile,daytype))
  save_table(bind_cols(tibble(term=rownames(kc)),as_tibble(kc),count_sample(d)[rep(1,nrow(kc)),]),paste0(label,'_phase_preliminary_k'),'results/diagnostics')
  save_table(patterns,paste0(label,'_phase_preliminary_patterns'),'results/diagnostics')
  estimated_rho<-pmin(.95,pmax(-.95,lag_cor(residuals(preliminary),d$AR_start)))
  rho<-if(is.null(rho_override))estimated_rho else rho_override
  rm(preliminary);gc()
  f<-fit_bam(d,formula,rho,paste(label,'phase AR1'))
  out<-list(fit=f,data=d,prep=prep,axis='phase',label=label,rho=rho,estimated_rho=estimated_rho,
    signature=signature,preliminary_k=kc)
  atomic_rds(out,path)
  out
}

daylight_support <- function(d,targets,x) {
  bind_rows(lapply(seq_len(nrow(targets)),function(i) {
    sub<-filter(d,abs(photo-targets$photo[i])<=1,quartile %in% c('Q1','Q4')) |>
      mutate(bin=floor((phase_display+12)/.5)) |>
      summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(bin,quartile,daytype))
    expand_grid(x=x,quartile=c('Q1','Q4'),daytype=c('Work','Free')) |>
      mutate(bin=floor((x+12)/.5)) |> left_join(sub,by=c('bin','quartile','daytype')) |>
      mutate(across(c(participants,participant_days,half_hours),~replace_na(.x,0L)),
        scenario=targets$scenario[i],photo_target=targets$photo[i])
  }))
}

daylight_matrices <- function(bundle,design) {
  x<-seq(-12,12-cfg$step,by=cfg$step); matrices<-list()
  for(i in seq_len(nrow(design$targets))) for(q in c('Q1','Q4')) {
    cell<-lapply(c('Work','Free'),function(day) {
      nd<-prediction_template(bundle,length(x));nd$t<-wrap24(x);nd$photo<-design$targets$photo[i]
      nd$quartile<-q;nd$direction<-if(q=='Q1')'Morning' else 'Evening';nd$extremity<-'Extreme';nd$daytype<-day
      # Full evaluation is essential: the new smooth has photo as its by variable.
      population_matrix(bundle,nd)
    })
    matrices[[paste(design$targets$scenario[i],q,sep='_')]]<-(1-design$free_weight)*cell[[1]]+design$free_weight*cell[[2]]
  }
  D<-lapply(design$targets$scenario,function(s)matrices[[paste0(s,'_Q4')]]-matrices[[paste0(s,'_Q1')]])
  names(D)<-design$targets$scenario;D$Longer_minus_shorter<-D$Longer-D$Shorter
  list(x=x,X=matrices,D=D)
}

daylight_basis_checks <- function(bundle) {
  f<-bundle$fit;ix<-which(vapply(f$smooth,function(sm)identical(sm$by,'photo')&&inherits(sm,'sz.interaction'),logical(1)))
  stopifnot(length(ix)==1L);sm<-f$smooth[[ix]];cols<-sm$first.para:sm$last.para
  grid<-expand_grid(t=c(0,6,12,18,24),quartile=levels(bundle$data$quartile))
  nd<-prediction_template(bundle,nrow(grid));nd$t<-grid$t;nd$quartile<-grid$quartile
  at<-function(photo) {z<-nd;z$photo<-photo;population_matrix(bundle,z)[,cols,drop=FALSE]}
  X0<-at(0);X1<-at(1);X2<-at(2)
  values<-c(length(cols),max(abs(rowsum(X2,grid$t))),max(abs(X0)),max(abs(X2-2*X1)),max(abs(X2[grid$t==0,]-X2[grid$t==24,])))
  checks<-tibble(check=c('Constrained cyclic sz coefficient count','Sum to zero across all quartiles','Zero moderation at centered photoperiod zero',
    'Linear scaling in centered photoperiod','Cyclic 0/24 endpoint agreement'),value=values,
    expected=c(33,0,0,0,0),tolerance=c(0,1e-12,1e-12,1e-12,1e-12),passed=abs(value-expected)<=tolerance)
  save_table(bind_cols(checks,count_sample(bundle$data)[rep(1,nrow(checks)),]),'daylight_basis_acceptance','results/provenance')
  stopifnot(all(checks$passed));checks
}

daylight_pointwise <- function(result,prefix='daylight') {
  z<-result$contrasts
  z<-z |> mutate(pointwise_lower=estimate-qnorm(.975)*se,
    pointwise_upper=estimate+qnorm(.975)*se,
    pointwise_ratio_lower=10^pointwise_lower,pointwise_ratio_upper=10^pointwise_upper,
    pointwise_excludes_one=pointwise_lower>0|pointwise_upper<0,
    pointwise_significant=supported&pointwise_excludes_one,
    pointwise_critical=qnorm(.975))
  result$contrasts<-z;result$reporting_version<-daylight_reporting_version
  result$information$pointwise_critical<-qnorm(.975)
  result$information$pointwise_interval<-'95% pointwise Gaussian; exploratory and unadjusted across the phase grid'
  save_table(z,paste0(prefix,'_contrasts'))
  save_table(z |> select(-starts_with('pointwise_')),paste0(prefix,'_simultaneous_contrasts'))
  save_table(z |> select(contrast,x,estimate,se,pointwise_lower,pointwise_upper,ratio,
    pointwise_ratio_lower,pointwise_ratio_upper,pointwise_excludes_one,pointwise_significant,
    supported,minimum_cell_participants,minimum_cell_days,morning_participants,morning_days,
    evening_participants,evening_days),paste0(prefix,'_pointwise_contrasts'))
  save_table(z |> filter(contrast=='Longer_minus_shorter') |> select(x,ratio,
    pointwise_ratio_lower,pointwise_ratio_upper,pointwise_excludes_one,pointwise_significant,
    supported,minimum_cell_participants,minimum_cell_days,morning_participants,morning_days,
    evening_participants,evening_days),paste0(prefix,'_pointwise_moderation_exclusions'))
  save_table(result$information,paste0(prefix,'_inference_information'))
  result
}

daylight_inference <- function(bundle,design) {
  prefix<-design$prefix
  path<-paste0('results/data/',prefix,'_inference.rds')
  signature<-digest::digest(list(bundle$signature,design$targets,design$free_weight,cfg$draws,cfg$step,daylight_version,daylight_support_n,daylight_reporting_version))
  if(file.exists(path)) {
    old<-readRDS(path)
    if(identical(old$signature,signature)) {
      changed<-!identical(old$reporting_version,daylight_reporting_version)
      old<-daylight_pointwise(old,prefix)
      if(changed)atomic_rds(old,path)
      return(old)
    }
    rm(old)
  }
  mm<-daylight_matrices(bundle,design);d<-bundle$data;f<-bundle$fit;V<-covariance_info(f);beta<-coef(f)
  support<-daylight_support(d,design$targets,mm$x)
  gate<-support |> summarise(minimum_cell_participants=min(participants),minimum_cell_days=min(participant_days),.by=c(scenario,x)) |>
    mutate(supported=minimum_cell_participants>=daylight_support_n)
  patterns<-bind_rows(lapply(names(mm$X),function(key) {
    parts<-strsplit(key,'_',fixed=TRUE)[[1]];sc<-parts[1];q<-parts[2]
    n<-count_sample(filter(d,quartile==q));local<-count_sample(filter(d,quartile==q,abs(photo-design$targets$photo[match(sc,design$targets$scenario)])<=1))
    bind_cols(tibble(scenario=sc,quartile=q,x=mm$x),linear_summary(mm$X[[key]],beta,V$V)) |>
      mutate(participants=n$participants,participant_days=n$participant_days,half_hours=n$half_hours,
        local_participants=local$participants,local_participant_days=local$participant_days,local_half_hours=local$half_hours) |>
      left_join(gate,by=c('scenario','x'))
  })) |> intensity_display()
  DD<-do.call(rbind,mm$D);keep<-which(colSums(abs(DD))>1e-12);ds<-linear_summary(DD,beta,V$V)
  set.seed(cfg$seed+701L)
  draws<-as.matrix(gratia::gaussian_draws(f,n=cfg$draws,unconditional=V$unconditional))
  centered<-sweep(draws[,keep,drop=FALSE],2,beta[keep],'-');maxz<-numeric(nrow(draws));usable<-ds$se>1e-10
  for(ii in split(seq_len(nrow(draws)),ceiling(seq_len(nrow(draws))/250))) {
    zz<-centered[ii,,drop=FALSE]%*%t(DD[usable,keep,drop=FALSE]);zz<-sweep(zz,2,ds$se[usable],'/');maxz[ii]<-apply(abs(zz),1,max)
  }
  critical<-unname(quantile(maxz,.95,type=8));rm(draws,centered);gc()
  contrasts<-bind_rows(lapply(names(mm$D),function(sc) {
    ss<-if(sc=='Longer_minus_shorter')gate |> filter(scenario %in% c('Shorter','Longer')) |>
      summarise(minimum_cell_participants=min(minimum_cell_participants),minimum_cell_days=min(minimum_cell_days),supported=all(supported),.by=x) else filter(gate,scenario==sc) |> select(-scenario)
    bind_cols(tibble(contrast=sc,x=mm$x),linear_summary(mm$D[[sc]],beta,V$V,critical)) |>
      left_join(ss,by='x') |> mutate(significant=supported&(lower>0|upper<0),critical=critical,
        morning_participants=n_distinct(d$participant_key[d$quartile=='Q1']),morning_days=n_distinct(d$participant_day_key[d$quartile=='Q1']),
        evening_participants=n_distinct(d$participant_key[d$quartile=='Q4']),evening_days=n_distinct(d$participant_day_key[d$quartile=='Q4']))
  })) |> contrast_display()
  info<-bind_cols(tibble(critical=critical,draws=cfg$draws,curves=length(mm$D),time_points=length(mm$x),coverage=.95,
    covariance=V$label,free_weight=design$free_weight),count_sample(d))
  save_table(patterns,paste0(prefix,'_patterns'));save_table(contrasts,paste0(prefix,'_contrasts'));save_table(support,paste0(prefix,'_local_support'))
  save_table(info,paste0(prefix,'_inference_information'))
  result<-daylight_pointwise(list(patterns=patterns,contrasts=contrasts,support=support,information=info,signature=signature),prefix)
  atomic_rds(result,path)
  result
}

daylight_windows <- function(bundle,inference,design,interval_type=c('pointwise','simultaneous')) {
  prefix<-design$prefix
  interval_type<-match.arg(interval_type)
  d<-bundle$data;out<-list();x<-sort(unique(inference$contrasts$x));n<-length(x)
  for(id in unique(inference$contrasts$contrast)) {
    z<-filter(inference$contrasts,contrast==id) |> arrange(x)
    if(interval_type=='pointwise')z$significant<-z$pointwise_significant
    sign<-ifelse(z$significant,ifelse(z$estimate>0,1L,-1L),0L)
    runs<-split(seq_len(n),cumsum(c(TRUE,diff(sign)!=0)))
    if(length(runs)>1&&sign[1]!=0&&sign[1]==sign[n]) {runs[[1]]<-c(runs[[length(runs)]],runs[[1]]);runs<-runs[-length(runs)]}
    targets<-if(id=='Longer_minus_shorter')c('Shorter','Longer')else id
    for(ix in runs) {
      if(!sign[ix[1]])next
      start<-wrap12(z$x[ix[1]]-cfg$step/2);duration<-length(ix)*cfg$step
      sub<-d |> filter(quartile %in% c('Q1','Q4'),wrap24(phase_display-start)<duration-1e-10)
      counts<-bind_rows(lapply(targets,function(sc) {
        sub |> filter(abs(photo-design$targets$photo[match(sc,design$targets$scenario)])<=1) |>
          summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=quartile) |>
          mutate(scenario=sc)
      }))
      if(nrow(counts)!=2L*length(targets)||any(counts$participants<daylight_support_n))next
      # Keep each scenario's denominator separate: days can overlap across targets.
      nr<-counts |> pivot_wider(names_from=c(scenario,quartile),values_from=c(participants,participant_days,half_hours))
      row<-tibble(contrast=id,interval_type=interval_type,direction=if(sign[ix[1]]>0)'Above one'else'Below one',start_phase=start,
        end_phase=wrap12(start+duration),duration_hours=duration,
        crosses_boundary=start>=wrap12(start+duration),minimum_ratio=min(z$ratio[ix]),maximum_ratio=max(z$ratio[ix]),
        geometric_mean_ratio=10^mean(z$estimate[ix]),minimum_cell_participants=min(z$minimum_cell_participants[ix]))
      out[[length(out)+1L]]<-bind_cols(row,nr)
    }
  }
  result<-if(length(out))bind_rows(out)else tibble(contrast=character(),direction=character(),start_phase=numeric(),end_phase=numeric(),duration_hours=numeric())
  save_table(result,if(interval_type=='pointwise')paste0(prefix,'_pointwise_windows')else paste0(prefix,'_significant_windows'))
  result
}

daylight_primary_intervals <- function(bundle,inference,design) {
  prefix<-design$prefix
  selected<-read.csv('results/tables/primary_phase_significant_windows.csv') |>
    filter(contrast=='evening_morning_Extreme_Overall')
  rows<-list();counts<-list();mm<-daylight_matrices(bundle,design)
  V<-covariance_info(bundle$fit);beta<-coef(bundle$fit)
  for(i in seq_len(nrow(selected))) {
    start<-selected$start[i];duration<-selected$duration_hours[i]
    for(id in unique(inference$contrasts$contrast)) {
      z<-filter(inference$contrasts,contrast==id,wrap24(x-start)<duration-1e-10)
      targets<-if(id=='Longer_minus_shorter')c('Shorter','Longer')else id
      nr<-bind_rows(lapply(targets,function(sc) {
        bundle$data |> filter(quartile %in% c('Q1','Q4'),wrap24(phase_display-start)<duration-1e-10,
          abs(photo-design$targets$photo[match(sc,design$targets$scenario)])<=1) |>
          summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=quartile) |>
          mutate(scenario=sc,interval=i,contrast=id)
      }))
      counts[[length(counts)+1L]]<-nr
      # Covariance of the averaged design vector includes correlation among
      # phases. The selected interval is treated as fixed; these are exploratory
      # Gaussian CIs without selection or multiplicity correction.
      ix<-which(wrap24(mm$x-start)<duration-1e-10)
      avg<-linear_summary(matrix(colMeans(mm$D[[id]][ix,,drop=FALSE]),nrow=1),beta,V$V)
      stopifnot(abs(avg$estimate-mean(z$estimate))<1e-9)
      rows[[length(rows)+1L]]<-tibble(interval=i,start_phase=start,end_phase=selected$end[i],duration_hours=duration,
        contrast=id,mean_log10_difference=avg$estimate,mean_log10_se=avg$se,ratio=10^avg$estimate,
        ratio_lower=10^avg$lower,ratio_upper=10^avg$upper,
        percent_change=100*(10^avg$estimate-1),percent_lower=100*(10^avg$lower-1),percent_upper=100*(10^avg$upper-1),
        simultaneous_envelope_lower=10^mean(z$lower),simultaneous_envelope_upper=10^mean(z$upper),
        interval_type='95% Gaussian CI for average log contrast; selected interval fixed; exploratory and unadjusted',
        covariance=V$label,grid_points=length(ix),supported=all(z$supported)&&all(nr$participants>=daylight_support_n),
        minimum_cell_participants=min(z$minimum_cell_participants),
        counts=paste(paste0(nr$scenario,' ',nr$quartile,': ',nr$participants,' participants / ',nr$participant_days,' days / ',nr$half_hours,' half-hours'),collapse='; '))
    }
  }
  out<-bind_rows(rows);save_table(out,paste0(prefix,'_primary_interval_summaries'))
  save_table(out |> select(interval,start_phase,end_phase,duration_hours,contrast,
    ratio,simultaneous_envelope_lower,simultaneous_envelope_upper,supported,minimum_cell_participants,counts),
    paste0(prefix,'_primary_interval_simultaneous_envelopes'))
  save_table(bind_rows(counts),paste0(prefix,'_primary_interval_samples'))
  out
}

daylight_plots <- function(bundle,inference,design) {
  prefix<-design$prefix
  d<-bundle$data;z<-inference$patterns
  pointwise_contrasts<-inference$contrasts |> mutate(ratio_lower=pointwise_ratio_lower,
    ratio_upper=pointwise_ratio_upper) |>
    group_by(contrast) |> arrange(x,.by_group=TRUE) |>
    mutate(sig_segment=cumsum(c(TRUE,diff(pointwise_significant)!=0))) |> ungroup()
  z$scenario<-factor(z$scenario,levels=c('Shorter','Typical','Longer'))
  counts<-z |> distinct(scenario,quartile,local_participants,local_participant_days) |>
    summarise(label=paste(paste0(ifelse(quartile=='Q1','Morning','Evening'),': ',local_participants,' participants / ',local_participant_days,' participant-days'),collapse='\n'),.by=scenario)
  labs<-setNames(paste0(LETTERS[1:3],'. ',design$targets$scenario,' days (',sprintf('%+.2f',design$targets$photo),' h)\n',counts$label[match(design$targets$scenario,counts$scenario)]),design$targets$scenario)
  totals<-d |> filter(quartile %in% c('Q1','Q4')) |> summarise(participants=n_distinct(participant_key),days=n_distinct(participant_day_key),.by=quartile)
  total_label<-paste(paste0(totals$quartile,': ',totals$participants,' people / ',totals$days,' days'),collapse='; ')
  p<-ggplot(z,aes(x,estimate_melEDI,color=quartile,fill=quartile))+
    primary_phase_interval_band()+dlmo_reference_band()+
    geom_ribbon(aes(ymin=lower_melEDI,ymax=upper_melEDI),alpha=.18,color=NA)+geom_line(linewidth=.7)+
    facet_wrap(~scenario,ncol=3,labeller=as_labeller(labs))+
    scale_y_continuous(transform=LightLogR::symlog_trans(base=10,thr=1,scale=1),breaks=c(0,1,10,100,1000))+
    scale_x_continuous(breaks=seq(-12,12,6),limits=c(-12,12))+
    scale_color_manual(values=c(Q1='#0072B2',Q4='#D55E00'),labels=c(Q1='Extreme morning',Q4='Extreme evening'))+
    scale_fill_manual(values=c(Q1='#0072B2',Q4='#D55E00'),labels=c(Q1='Extreme morning',Q4='Extreme evening'))+
    labs(title='Exploratory circadian phase curves across photoperiod',subtitle=sample_label(d),x='Hours relative to MSFsc',y='Melanopic EDI (lx; symlog axis)',color=NULL,fill=NULL,
      caption=paste0('Targets are centered within site. Panel counts: observed days within ±1 h of each target. Pointwise 95% bands.\n',
        total_label,'. Grey windows with dashed boundaries: primary phase intervals. Striped band: DLMO reference.\n',
        'Curves average the same work/free-day mix and exclude site, participant and day effects. Bands remain exploratory.'))
  save_plot(p,paste0(prefix,'_phase_patterns'),14,6)
  zz<-pointwise_contrasts |> filter(contrast!='Longer_minus_shorter') |> mutate(scenario=factor(contrast,levels=levels(z$scenario)))
  ratio_limits<-range(c(1,zz$ratio,zz$ratio_lower[zz$supported],zz$ratio_upper[zz$supported]))
  # Retain unsupported estimates, but show inferential bands only with local support.
  p<-ggplot(zz,aes(x,ratio))+
    primary_phase_interval_band()+dlmo_reference_band()+
    geom_hline(yintercept=1,color='grey50')+
    geom_ribbon(aes(ymin=ifelse(supported,ratio_lower,NA),ymax=ifelse(supported,ratio_upper,NA)),fill='#6A3D9A',alpha=.18)+
    geom_line(color='#6A3D9A',linewidth=.7)+
    geom_line(data=filter(zz,pointwise_significant),aes(group=interaction(contrast,sig_segment)),color='#6A3D9A',linewidth=1.6)+geom_point(data=filter(zz,!supported),shape=4,color='grey50',size=.8)+
    facet_wrap(~scenario,ncol=3,labeller=as_labeller(labs))+
    scale_y_log10(limits=ratio_limits,oob=scales::oob_keep,breaks=c(.01,.03,.1,.3,1,3,10,30,100),labels=function(v)paste0(v,'×'))+
    scale_x_continuous(breaks=seq(-12,12,6),limits=c(-12,12))+
    labs(title='Exploratory extreme-evening / extreme-morning ratio by circadian phase',subtitle=sample_label(d),x='Hours relative to MSFsc',y='Ratio of melEDI + 0.1',
      caption='95% pointwise bands; exploratory comparisons without adjustment across phases or scenarios.\nBands are withheld where any extreme-group × day-type cell has fewer than four people in the local phase bin and ±1 h of the target.\nBold segments: supported 95% intervals exclude 1. Crosses mark unsupported estimates. Panel counts are local covariate counts.\nGrey windows with dashed boundaries: primary phase intervals. Striped band: DLMO reference.')
  save_plot(p,paste0(prefix,'_phase_contrasts'),14,6)
  zz<-filter(pointwise_contrasts,contrast=='Longer_minus_shorter')
  ratio_limits<-range(c(1,zz$ratio,zz$ratio_lower[zz$supported],zz$ratio_upper[zz$supported]))
  p<-ggplot(zz,aes(x,ratio))+primary_phase_interval_band()+dlmo_reference_band()+
    geom_hline(yintercept=1,color='grey50')+geom_ribbon(aes(ymin=ifelse(supported,ratio_lower,NA),ymax=ifelse(supported,ratio_upper,NA)),fill='#6A3D9A',alpha=.18)+
    geom_line(color='#6A3D9A',linewidth=.7)+
    geom_line(data=filter(zz,pointwise_significant),aes(group=interaction(contrast,sig_segment)),color='#6A3D9A',linewidth=1.6)+geom_point(data=filter(zz,!supported),shape=4,color='grey50',size=.8)+
    scale_y_log10(limits=ratio_limits,oob=scales::oob_keep,breaks=c(.01,.03,.1,.3,1,3,10,30,100),labels=function(v)paste0(v,'×'))+
    scale_x_continuous(breaks=seq(-12,12,3),limits=c(-12,12))+
    labs(title='Does the extreme-group ratio change with daylength?',subtitle=sample_label(d),x='Hours relative to MSFsc',y='Longer / shorter ratio of ratios',
      caption=paste0('Above 1: the evening/morning ratio is larger on longer days; below 1: it is smaller. Circadian phase is relative to MSFsc.\n',
        'Shorter: ',gsub('\n','; ',counts$label[counts$scenario=='Shorter']),'.\n',
        'Longer: ',gsub('\n','; ',counts$label[counts$scenario=='Longer']),'.\n',
        '95% pointwise exploratory band, unadjusted across phase. Bold segments: supported 95% intervals exclude 1. Crosses mark unsupported phase bins.\n',
        'Grey windows with dashed boundaries: primary phase intervals. Striped band: DLMO reference.'))
  save_plot(p,paste0(prefix,'_phase_moderation'),11,6)
}

daylight_solar <- function(bundle) {
  d<-bundle$data;context<-readRDS(input_paths()[['glasses_context']])
  keys<-c('site','Id','position','local_date','clock_bin')
  solar<-context |> select(all_of(keys),civil_dawn_wall_minute,civil_dusk_wall_minute,solar_depression_deg)
  stopifnot(!anyDuplicated(solar[keys]))
  joined<-left_join(d,solar,by=keys)
  stopifnot(nrow(joined)==nrow(d),!anyNA(joined$civil_dawn_wall_minute),!anyNA(joined$civil_dusk_wall_minute))
  days<-joined |> distinct(participant_key,participant_day_key,site,quartile,msf_hours,photoperiod_hours,
    civil_dawn_wall_minute,civil_dusk_wall_minute,solar_depression_deg)
  stopifnot(!anyDuplicated(days$participant_day_key),all(days$solar_depression_deg==6))
  reference<-days |> group_by(participant_key) |> mutate(w=1/n()) |> ungroup()
  median<-daylight_weighted_quantile(reference$photoperiod_hours,reference$w,.5)
  days<-days |> mutate(daylength=if_else(photoperiod_hours<=median,'Shorter','Longer'))
  extreme<-filter(days,quartile %in% c('Q1','Q4')) |>
    group_by(quartile,daylength,participant_key) |> mutate(w=1/n()) |> ungroup() |>
    group_by(quartile,daylength) |> mutate(w=w/n_distinct(participant_key)) |> ungroup()
  x<-seq(-12,12-cfg$step,by=cfg$step)
  estimates<-bind_rows(lapply(x,function(xx) {
    extreme |> mutate(clock_minute=wrap24(xx+msf_hours)*60,
      civil_available=((clock_minute-civil_dawn_wall_minute)%%1440)<=((civil_dusk_wall_minute-civil_dawn_wall_minute)%%1440)) |>
      summarise(availability=sum(w*civil_available),participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),.by=c(quartile,daylength)) |> mutate(x=xx)
  }))
  counts<-joined |> left_join(select(days,participant_day_key,daylength),by='participant_day_key') |>
    filter(quartile %in% c('Q1','Q4')) |> summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),
      minimum_daylength=min(photoperiod_hours),maximum_daylength=max(photoperiod_hours),.by=c(quartile,daylength))
  stopifnot(all(estimates$availability>=-1e-12),all(estimates$availability<=1+1e-12))
  selected<-read.csv('results/tables/primary_phase_significant_windows.csv') |>
    filter(contrast=='evening_morning_Extreme_Overall') |> arrange(start) |>
    transmute(interval=row_number(),interval_label=if_else(end<0,'Before DLMO reference','After MSFsc'),
      start_phase=start,end_phase=end,duration_hours)
  stopifnot(nrow(selected)==2L)
  observed<-joined |> left_join(select(days,participant_day_key,daylength),by='participant_day_key') |>
    filter(quartile %in% c('Q1','Q4'))
  reference_summary<-bind_rows(lapply(seq_len(nrow(selected)),function(i) {
    interval<-selected[i,];ix<-which(wrap24(x-interval$start_phase)<interval$duration_hours-1e-10)
    observed_counts<-observed |> filter(wrap24(phase_display-interval$start_phase)<interval$duration_hours-1e-10) |>
      summarise(interval_observed_participants=n_distinct(participant_key),
        interval_observed_days=n_distinct(participant_day_key),interval_observed_half_hours=n(),.by=c(quartile,daylength))
    estimates |> filter(x %in% .env$x[ix]) |>
      summarise(mean_availability=mean(availability),minimum_availability=min(availability),
        maximum_availability=max(availability),.by=c(quartile,daylength)) |>
      left_join(counts,by=c('quartile','daylength')) |>
      left_join(observed_counts,by=c('quartile','daylength')) |>
      mutate(interval=.env$interval$interval,interval_label=.env$interval$interval_label,
        start_phase=.env$interval$start_phase,end_phase=.env$interval$end_phase,duration_hours=.env$interval$duration_hours,
        grid_points=length(ix),grid_day_evaluations=participant_days*length(ix),
        mean_percent=100*mean_availability,minimum_percent=100*minimum_availability,maximum_percent=100*maximum_availability,
        denominator='All retained days in this quartile/absolute-daylength cell with solar metadata; equal participant and within-participant day weights')
  }))
  save_table(estimates,'daylight_civil_availability');save_table(counts,'daylight_civil_samples')
  save_table(reference_summary,'daylight_civil_reference_summary')
  save_table(reference_summary,'daylight_civil_primary_intervals')
  save_table(selected,'daylight_civil_interval_definitions')
  save_table(bind_cols(tibble(weighted_median_hours=median,solar_depression_degrees=6,definition='Civil dawn-to-dusk; astronomical availability, not measured light source'),count_sample(d)),'daylight_civil_definition')
  labs<-counts |> arrange(quartile) |> summarise(label=paste(paste0(ifelse(quartile=='Q1','Morning','Evening'),': ',participants,' participants / ',participant_days,' participant-days'),collapse='\n'),.by=daylength)
  labels<-setNames(paste0(ifelse(labs$daylength=='Shorter','A. ','B. '),labs$daylength,' absolute days\n',labs$label),labs$daylength)
  estimates$daylength<-factor(estimates$daylength,levels=c('Shorter','Longer'))
  p<-ggplot(estimates,aes(x,availability,color=quartile))+
    primary_phase_interval_band()+dlmo_reference_band()+
    geom_line(linewidth=.8)+
    facet_wrap(~daylength,labeller=as_labeller(labels))+
    scale_y_continuous(labels=scales::label_percent(),limits=c(0,1),oob=scales::oob_keep)+scale_x_continuous(breaks=seq(-12,12,3),limits=c(-12,12))+
    scale_color_manual(values=c(Q1='#0072B2',Q4='#D55E00'),labels=c(Q1='Extreme morning',Q4='Extreme evening'))+
    labs(title='Civil daylight availability across the two primary phase windows',subtitle=sample_label(filter(d,quartile %in% c('Q1','Q4'))),x='Hours relative to MSFsc',y='Weighted share of days with civil daylight',color=NULL,
      caption=sprintf("Absolute daylength is split at %.2f h, using the full primary sample (%s).\nWithin each panel and group, participants receive equal weight and their days receive equal weight.\nGrey windows with dashed boundaries: before DLMO (%.2f to %.2f h) and after MSFsc (%.2f to %.2f h). Striped band: DLMO reference.\nCivil dawn-to-dusk uses solar altitude -6°. Counts include all days with solar metadata in each group, regardless of observed light at that phase.",
        median,sample_label(d),selected$start_phase[1],selected$end_phase[1],selected$start_phase[2],selected$end_phase[2]))
  save_plot(p,'daylight_civil_availability',12,6)
  list(estimates=estimates,counts=counts,median=median,reference_summary=reference_summary)
}

daylight_validation <- function(prepared) {
  raw<-prepared$frames$glasses;registry<-read.csv('results/validation/participant_folds.csv')
  stopifnot(!anyDuplicated(registry$participant_key),setequal(registry$participant_key,raw$participant_key))
  predictions<-samples<-list()
  for(i in 1:5) {
    id<-paste0('new_participant_fold_',i);assessment<-raw$participant_key %in% registry$participant_key[registry$fold==i]
    train<-raw[!assessment,];test<-raw[assessment,];prep<-fit_preprocess(train)
    stopifnot(!length(intersect(train$participant_key,test$participant_key)),!length(intersect(train$participant_day_key,test$participant_day_key)))
    signature<-digest::digest(list(train,test,prep,daylight_version,deparse(daylight_formula())))
    path<-paste0('results/validation/daylight_',id,'_phase.rds');out<-if(file.exists(path))readRDS(path)else NULL
    base<-readRDS(paste0('results/validation/',id,'_phase.rds'))
    stopifnot(isTRUE(all.equal(prep$cutoffs,base$cutoffs)),isTRUE(all.equal(prep$centers,base$centers)))
    if(is.null(out)||!identical(out$signature,signature)) {
      b<-daylight_fit(train,prep,label=paste0('daylight_',id))
      pr<-predict_assessment(b,test,'New participants') |> mutate(target='New participants',split=id,model='Photoperiod moderation')
      out<-list(predictions=pr,signature=signature,training=count_sample(train),assessment=count_sample(test),cutoffs=prep$cutoffs,centers=prep$centers,rho=b$rho,diagnostics=model_summary(b))
      atomic_rds(out,path);rm(b);gc()
    }
    baseline<-base$predictions |> mutate(model='Primary phase')
    key<-c('participant_key','participant_day_key','site','local_date','clock_bin')
    a<-arrange(baseline,across(all_of(key)));b<-arrange(out$predictions,across(all_of(key)))
    stopifnot(identical(a[key],b[key]),isTRUE(all.equal(a$observed,b$observed)),isTRUE(all.equal(a$baseline,b$baseline)))
    predictions[[id]]<-bind_rows(baseline,out$predictions)
    samples[[id]]<-bind_rows(lapply(list(base,out),function(o) tibble(split=id,training_participants=o$training$participants,training_days=o$training$participant_days,
      training_half_hours=o$training$half_hours,assessment_participants=o$assessment$participants,assessment_days=o$assessment$participant_days,
      assessment_half_hours=o$assessment$half_hours,cutoff25=o$cutoffs[1],cutoff50=o$cutoffs[2],cutoff75=o$cutoffs[3],rho=o$rho))) |>
      mutate(model=c('Primary phase','Photoperiod moderation'))
  }
  predictions<-bind_rows(predictions);samples<-bind_rows(samples)
  scores<-predictions |> group_by(model) |> group_modify(~score_rows(.x)) |> ungroup()
  folds<-predictions |> group_by(split,model) |> group_modify(~score_rows(.x)) |> ungroup()
  pp<-predictions |> summarise(mse=mean((predicted-observed)^2),mae=mean(abs(predicted-observed)),bias=mean(predicted-observed),null_mse=mean((baseline-observed)^2),.by=c(model,site,participant_key))
  a<-filter(pp,model=='Primary phase') |> arrange(site,participant_key);b<-filter(pp,model=='Photoperiod moderation') |> arrange(site,participant_key)
  stopifnot(identical(a$participant_key,b$participant_key),isTRUE(all.equal(a$null_mse,b$null_mse)))
  score<-function(z,ix)c(RMSE=sqrt(mean(z$mse[ix])),MAE=mean(z$mae[ix]),bias=mean(z$bias[ix]),predictive_R2=1-mean(z$mse[ix])/mean(z$null_mse[ix]))
  delta<-score(b,seq_len(nrow(b)))-score(a,seq_len(nrow(a)));strata<-split(seq_len(nrow(a)),a$site)
  set.seed(cfg$seed+702L)
  draws<-replicate(cfg$bootstrap,{ix<-unlist(lapply(strata,function(j)j[sample.int(length(j),length(j),replace=TRUE)]));score(b,ix)-score(a,ix)})
  ci<-apply(draws,1,quantile,c(.025,.975),type=8)
  differences<-bind_cols(tibble(metric=names(delta),expanded_minus_primary=as.numeric(delta),lower=ci[1,],upper=ci[2,],bootstrap_replicates=cfg$bootstrap),count_sample(raw)[rep(1,4),])
  save_table(predictions,'daylight_held_out_predictions','results/validation');save_table(samples,'daylight_validation_samples')
  save_table(scores,'daylight_validation_scores');save_table(folds,'daylight_validation_fold_scores');save_table(differences,'daylight_validation_paired_differences')
  list(scores=scores,samples=samples,differences=differences)
}

run_daylight <- function(prepared,validation=TRUE) {
  base_files<-c('results/models/primary_phase.rds','results/models/primary_clock.rds',
    'results/data/primary_phase_inference.rds','results/data/primary_clock_inference.rds',
    'results/validation/participant_folds.csv','results/validation/held_out_predictions.csv')
  before<-tools::md5sum(base_files)
  raw<-prepared$frames$glasses;prep<-fit_preprocess(raw,prepared$cutoffs);d<-apply_preprocess(raw,prep,'phase')
  design<-daylight_design(d);b<-daylight_fit(raw,prep);basis_checks<-daylight_basis_checks(b);inference<-daylight_inference(b,design)
  diagnostic_path<-'results/diagnostics/daylight_bundle.rds'
  diag<-if(file.exists(diagnostic_path))readRDS(diagnostic_path)else NULL
  if(is.null(diag)||!identical(diag$signature,b$signature)) {
    diag<-list(signature=b$signature,diagnostics=diagnose_model(b));atomic_rds(diag,diagnostic_path)
  }
  save_table(model_summary(b),'daylight_model_summary')
  windows<-daylight_windows(b,inference,design)
  simultaneous_windows<-daylight_windows(b,inference,design,interval_type='simultaneous')
  selected_intervals<-daylight_primary_intervals(b,inference,design)
  daylight_plots(b,inference,design);sensitivities<-daylight_sensitivities(b);solar<-daylight_solar(b)
  checks<-tibble(check=c('Identical primary analysis observations','Preserved true-elapsed AR ordering','All target-level quartile/day cells have four participants',
    'All significant grid points pass local support','Scenario contrast signs and longer-minus-shorter identity','No site has five participants in both extreme quartiles'),
    passed=c(identical(d,b$data),identical(d$AR_start,b$data$AR_start)&&identical(d$source_utc_start,b$data$source_utc_start),
      all(design$support$participants>=daylight_support_n),all(!inference$contrasts$significant|inference$contrasts$supported),
      {z<-inference$contrasts;max(abs(filter(z,contrast=='Longer_minus_shorter')$estimate-(filter(z,contrast=='Longer')$estimate-filter(z,contrast=='Shorter')$estimate)))<1e-9},
      {z<-design$site |> summarise(overlap=all(participants>=5),.by=site);!any(z$overlap)}))
  rm(b);gc()
  val<-if(validation)daylight_validation(prepared)else NULL
  after<-tools::md5sum(base_files)
  checks<-bind_rows(checks,select(basis_checks,check,passed),tibble(check='Primary model, inference and validation files unchanged',passed=identical(before,after)))
  checks<-bind_rows(checks,tibble(check=c('Pointwise moderation covers all 144 phase grid points',
    'Pointwise exclusions obey local support','Civil summaries use the two primary phase windows'),
    passed=c(sum(inference$contrasts$contrast=='Longer_minus_shorter')==144L,
      all(!inference$contrasts$pointwise_significant|inference$contrasts$supported),
      n_distinct(solar$reference_summary$interval)==2L&&nrow(solar$reference_summary)==8L)))
  save_table(checks,'daylight_acceptance_checks','results/provenance')
  save_table(tibble(file=names(before),before_md5=unname(before),after_md5=unname(after),unchanged=before==after),'daylight_primary_cache_preservation','results/provenance')
  source_files<-paste0('R/',c('setup','data','models','inference','diagnostics','validation','daylight'),'.R')
  save_table(tibble(file=source_files,md5=unname(tools::md5sum(source_files))),
    'daylight_source_checksums','results/provenance')
  manifest<-input_manifest()
  save_table(manifest |> select(role,snapshot_path,upstream_path,source_repository,source_commit,md5,sha256,bytes),
    'daylight_input_provenance','results/provenance')
  writeLines(c(paste('R:',getRversion()),paste('Daylight implementation:',daylight_version),
    paste('Daylight reporting:',daylight_reporting_version),
    paste('Frozen preprocessed observation digest:',digest::digest(d)),
    'Command: Rscript --vanilla run.R daylight',
    'Model fitting: full phase model plus constrained cyclic quartile-specific linear photoperiod change.',
    'Five participant-held-out splits read from the unchanged primary registry; training-only preprocessing and AR estimation.',
    'All input data are local snapshots described in daylight_input_provenance.csv.',
    paste('Inference covariance:',inference$information$covariance)),
    'results/provenance/daylight_run.txt')
  capture.output(sessionInfo(),file='results/provenance/daylight_session.txt')
  stopifnot(all(checks$passed))
  invisible(list(design=design,inference=inference,windows=windows,simultaneous_windows=simultaneous_windows,
    selected_intervals=selected_intervals,sensitivities=sensitivities,solar=solar,diagnostics=diag,validation=val))
}

# Reuse the fitted moderation model for both target sensitivities. These change
# the prediction contrast, not the model, sample, site centers or day-type mix.
daylight_sensitivities <- function(bundle) {
  out<-list()
  for(specification in c('quartiles','p10_p90')) {
    design<-daylight_design(bundle$data,specification,paste0('daylight_',specification))
    inference<-daylight_inference(bundle,design)
    intervals<-daylight_primary_intervals(bundle,inference,design)
    daylight_windows(bundle,inference,design)
    daylight_plots(bundle,inference,design)
    out[[specification]]<-list(design=design,inference=inference,intervals=intervals)
  }
  specs<-c('fixed','quartiles','p10_p90')
  prefixes<-c('daylight','daylight_quartiles','daylight_p10_p90')
  targets<-bind_rows(lapply(seq_along(specs),function(i)
    read.csv(paste0('results/tables/',prefixes[i],'_scenario_targets.csv')) |> mutate(specification=specs[i])))
  intervals<-bind_rows(lapply(seq_along(specs),function(i)
    read.csv(paste0('results/tables/',prefixes[i],'_primary_interval_summaries.csv')) |>
      filter(contrast=='Longer_minus_shorter') |> mutate(specification=specs[i])))
  save_table(targets,'daylight_sensitivity_targets');save_table(intervals,'daylight_sensitivity_intervals')
  out
}
