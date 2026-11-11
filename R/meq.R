# MEQ sensitivity: pooled morning/evening groups and a same-sample MSFsc comparator.
# The existing model response, cyclic bases, photoperiod standardization and AR resets are retained.
# Run from the project root: Rscript --vanilla R/meq.R

meq_prepare <- function(prepared) {
  d<-prepared$frames$glasses
  category_order<-c('Definitely morning type','Moderately morning type','Intermediate',
                    'Moderately evening type','Definitely evening type')
  expected<-case_when(d$meq>=70~category_order[1],d$meq>=59~category_order[2],
    d$meq>=42~category_order[3],d$meq>=31~category_order[4],d$meq>=16~category_order[5],TRUE~NA_character_)
  stopifnot(all(d$meq_type==expected),all(d$meq>=16 & d$meq<=86))
  categories<-d |> summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),
    half_hours=n(),minimum_score=min(meq),maximum_score=max(meq),.by=meq_type) |>
    mutate(meq_type=factor(meq_type,levels=category_order)) |> arrange(meq_type)
  save_table(categories,'meq_category_sample')
  save_table(d |> summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),
    half_hours=n(),.by=c(meq_type,daytype)),'meq_category_day_sample')
  d<-d |> mutate(meq_direction=case_when(grepl('morning',meq_type,fixed=TRUE)~'Morning',
    grepl('evening',meq_type,fixed=TRUE)~'Evening',TRUE~NA_character_),
    original_quartile=as.character(quartile),msf_direction=as.character(direction))
  raw<-filter(d,!is.na(meq_direction))
  flow<-bind_rows(mutate(count_sample(d),stage='Complete near-eye sample'),
    mutate(count_sample(filter(d,!is.na(meq))),stage='MEQ available'),
    mutate(count_sample(filter(d,meq_type=='Intermediate')),stage='Excluded: intermediate MEQ'),
    mutate(count_sample(raw),stage='MEQ morning/evening sensitivity'))
  save_table(flow,'meq_sample_flow')
  save_table(raw |> select(participant_key,participant_day_key,site,local_date,clock_bin,clock,phase_display,
    meq,meq_type,meq_direction,msf_direction,original_quartile,daytype,msf_hours,photoperiod_hours,
    AR_start,true_elapsed_sequence_id),'meq_sample_registry')
  save_table(raw |> summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),
    half_hours=n(),.by=c(meq_direction,msf_direction)),
    'meq_classification_overlap')
  site_support<-raw |> summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),
    half_hours=n(),.by=c(site,meq_direction))
  save_table(site_support,'meq_site_support')
  raw
}

meq_fit <- function(raw,prep,axis,classification,k=12L,kp=10L) {
  d<-apply_preprocess(raw,prep,axis)
  group<-if(classification=='MEQ')d$meq_direction else d$msf_direction
  d$quartile<-factor(group,levels=c('Morning','Evening'))
  d$direction<-d$quartile
  label<-if(classification=='MEQ')'meq_defined'else'meq_matched_msfsc'
  formula<-model_formula(k=k,kp=kp)
  signature<-digest::digest(list(data=d,formula=deparse(formula),classification=classification,version='meq-1.0.0',mgcv=as.character(packageVersion('mgcv'))))
  path<-paste0('results/models/',label,'_',axis,'.rds')
  if(file.exists(path)) {old<-readRDS(path);if(identical(old$signature,signature))return(old);rm(old);gc()}
  pre<-fit_bam(d,formula,0,paste(label,axis,'preliminary'))
  set.seed(cfg$seed+340L)
  kc<-mgcv::k.check(pre,subsample=5000,n.rep=100)
  save_table(bind_cols(tibble(term=rownames(kc)),as_tibble(kc),count_sample(d)[rep(1,nrow(kc)),]),
    paste0(label,'_',axis,'_preliminary_k'),'results/diagnostics')
  save_table(d |> mutate(residual=residuals(pre),hour=floor(t)) |>
    summarise(mean_residual=mean(residual),sd_residual=sd(residual),participants=n_distinct(participant_key),
      participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(quartile,daytype,hour)),
    paste0(label,'_',axis,'_preliminary_patterns'),'results/diagnostics')
  rho<-pmin(.95,pmax(-.95,lag_cor(residuals(pre),d$AR_start)))
  rm(pre);gc()
  f<-fit_bam(d,formula,rho,paste(label,axis,'AR1'))
  out<-list(fit=f,data=d,prep=prep,axis=axis,label=label,classification=classification,rho=rho,
    signature=signature,preliminary_k=kc)
  atomic_rds(out,path)
  out
}

meq_cells <- function(bundle) {
  ref<-reference_days(bundle$data) |> summarise(w=sum(w),.by=photo)
  x<-if(bundle$axis=='clock')seq(0,24-cfg$step,by=cfg$step)else seq(-12,12-cfg$step,by=cfg$step)
  nd<-prediction_template(bundle,length(x));nd$t<-wrap24(x);nd$photo<-0
  base<-population_matrix(bundle,nd)
  for(sm in bundle$fit$smooth)if('photo'%in%sm$term) {
    j<-sm$first.para:sm$last.para;avg<-matrix(0,length(x),length(j))
    for(i in seq_len(nrow(ref))) {nd$photo<-ref$photo[i];avg<-avg+ref$w[i]*mgcv::PredictMat(sm,nd)}
    base[,j]<-avg
  }
  nd$photo<-0;delta<-base-population_matrix(bundle,nd)
  cells<-expand_grid(group=c('Morning','Evening'),daytype=c('Work','Free'))
  X<-lapply(seq_len(nrow(cells)),function(i) {
    z<-nd;z$quartile<-cells$group[i];z$direction<-cells$group[i];z$daytype<-cells$daytype[i]
    population_matrix(bundle,z)+delta
  })
  names(X)<-paste(cells$group,cells$daytype,sep='_')
  list(x=x,X=X)
}

meq_support <- function(d,group,context,axis,x) {
  z<-filter(d,as.character(quartile)==.env$group)
  if(context!='Overall')z<-filter(z,as.character(daytype)==.env$context)
  origin<-if(axis=='clock')0 else -12
  local<-z |> mutate(bin=floor(((if(axis=='clock')clock else phase_display)-origin)/.5)) |>
    summarise(n=n_distinct(participant_key),days=n_distinct(participant_day_key),.by=bin)
  j<-match(floor((x-origin)/.5),local$bin)
  list(total=count_sample(z),n=replace_na(local$n[j],0L),days=replace_na(local$days[j],0L))
}

meq_inference <- function(bundle) {
  d<-bundle$data;cells<-meq_cells(bundle);x<-cells$x;V<-covariance_info(bundle$fit);beta<-coef(bundle$fit)
  ref<-reference_days(d);free_weight<-sum(ref$w[ref$daytype=='Free'])
  patterns<-D<-list();curves<-list()
  for(context in c('Overall','Work','Free'))for(group in c('Morning','Evening')) {
    key<-paste(group,context,sep='_')
    X<-if(context=='Overall')(1-free_weight)*cells$X[[paste0(group,'_Work')]]+free_weight*cells$X[[paste0(group,'_Free')]]else cells$X[[key]]
    patterns[[key]]<-X;su<-meq_support(d,group,context,bundle$axis,x)
    curves[[key]]<-bind_cols(tibble(classification=bundle$classification,axis=bundle$axis,x=x,context=context,group=group),linear_summary(X,beta,V$V)) |>
      mutate(participants=su$total$participants,participant_days=su$total$participant_days,half_hours=su$total$half_hours,
        support_participants=su$n,support_days=su$days)
  }
  for(context in c('Overall','Work','Free'))D[[context]]<-patterns[[paste0('Evening_',context)]]-patterns[[paste0('Morning_',context)]]
  DD<-do.call(rbind,D);s<-linear_summary(DD,beta,V$V);keep<-which(colSums(abs(DD))>1e-12)
  set.seed(cfg$seed+350L+as.integer(bundle$axis=='phase')+2L*as.integer(bundle$classification!='MEQ'))
  draws<-as.matrix(gratia::gaussian_draws(bundle$fit,n=cfg$draws,unconditional=V$unconditional))
  centered<-sweep(draws[,keep,drop=FALSE],2,beta[keep],'-');maxz<-numeric(nrow(draws));usable<-s$se>1e-10
  for(ii in split(seq_len(nrow(draws)),ceiling(seq_len(nrow(draws))/250))) {
    zz<-centered[ii,,drop=FALSE]%*%t(DD[usable,keep,drop=FALSE]);zz<-sweep(zz,2,s$se[usable],'/');maxz[ii]<-apply(abs(zz),1,max)
  }
  # Apply the same within-coordinate convention separately to each
  # classification sensitivity, covering all three day-type contexts.
  critical<-unname(quantile(maxz,cfg$curve_coverage,type=8))
  contrasts<-list()
  for(context in names(D)) {
    a<-meq_support(d,'Evening',context,bundle$axis,x);b<-meq_support(d,'Morning',context,bundle$axis,x)
    contexts<-if(context=='Overall')c('Work','Free')else context
    local<-lapply(contexts,function(ct)pmin(meq_support(d,'Morning',ct,bundle$axis,x)$n,meq_support(d,'Evening',ct,bundle$axis,x)$n))
    support<-Reduce(pmin,local)
    contrasts[[context]]<-bind_cols(tibble(classification=bundle$classification,axis=bundle$axis,x=x,context=context),linear_summary(D[[context]],beta,V$V,critical)) |>
      mutate(evening_participants=a$total$participants,evening_days=a$total$participant_days,evening_half_hours=a$total$half_hours,
        morning_participants=b$total$participants,morning_days=b$total$participant_days,morning_half_hours=b$total$half_hours,
        minimum_cell_participants=support,supported=support>=cfg$support_n,
        candidate_interval=supported & (lower>0|upper<0)) |> contrast_display()
  }
  stem<-paste(bundle$label,bundle$axis,sep='_')
  out<-list(curves=intensity_display(bind_rows(curves)),contrasts=bind_rows(contrasts),critical=critical,covariance=V$label,free_weight=free_weight)
  out$windows<-meq_windows(out$contrasts,d,bundle$axis)
  save_table(out$curves,paste0(stem,'_curves'));save_table(out$contrasts,paste0(stem,'_contrasts'))
  save_table(out$windows,paste0(stem,'_candidate_windows'))
  save_table(bind_cols(tibble(classification=bundle$classification,axis=bundle$axis,critical=critical,draws=cfg$draws,
    covariance=V$label,within_fit_coverage=cfg$curve_coverage,
    inference_family='All three day contexts and grid points within this classification and coordinate',
    across_coordinate_adjustment='None; separate time-coordinate and classification sensitivity questions',
    minimum_cell_support=cfg$support_n,free_weight=free_weight),count_sample(d)),paste0(stem,'_inference_metadata'))
  out
}

meq_windows <- function(contrasts,d,axis) {
  out<-list();origin<-if(axis=='clock')0 else -12
  for(context in unique(contrasts$context)) {
    z<-filter(contrasts,.data$context==.env$context) |> arrange(x)
    sg<-ifelse(z$candidate_interval,ifelse(z$estimate>0,1L,-1L),0L)
    rr<-rle(sg);ends<-cumsum(rr$lengths);groups<-Map(seq.int,ends-rr$lengths+1L,ends);values<-rr$values
    if(length(groups)>1 && values[1]!=0 && values[1]==tail(values,1)) {
      groups[[1]]<-c(tail(groups,1)[[1]],groups[[1]]);groups<-head(groups,-1);values<-head(values,-1)
    }
    for(j in which(values!=0)) {
      ix<-groups[[j]];start<-wrap24(z$x[ix[1]]-cfg$step/2-origin)+origin;duration<-length(ix)*cfg$step;end<-wrap24(start+duration-origin)+origin
      coordinate<-if(axis=='clock')d$clock else d$phase_display
      rows<-d[wrap24(coordinate-start)<duration-1e-10,]
      if(context!='Overall')rows<-filter(rows,as.character(daytype)==.env$context)
      a<-count_sample(filter(rows,as.character(quartile)=='Evening'));b<-count_sample(filter(rows,as.character(quartile)=='Morning'))
      if(min(a$participants,b$participants)<cfg$support_n)next
      out[[length(out)+1L]]<-tibble(classification=z$classification[1],axis=axis,context=context,start=start,end=end,
        wraps_boundary=start>=end,duration_hours=duration,direction=if(values[j]>0)'Evening higher'else'Evening lower',
        mean_log10_difference=mean(z$estimate[ix]),factor_at_mean_difference=10^mean(z$estimate[ix]),
        evening_participants=a$participants,evening_days=a$participant_days,evening_half_hours=a$half_hours,
        morning_participants=b$participants,morning_days=b$participant_days,morning_half_hours=b$half_hours)
    }
  }
  if(length(out))bind_rows(out)else tibble(classification=character(),axis=character(),context=character(),start=double(),end=double(),
    wraps_boundary=logical(),duration_hours=double(),direction=character(),mean_log10_difference=double(),factor_at_mean_difference=double(),
    evening_participants=integer(),evening_days=integer(),evening_half_hours=integer(),morning_participants=integer(),morning_days=integer(),morning_half_hours=integer())
}

meq_pattern_figure <- function(curves,classification) {
  panels<-list();zall<-filter(curves,.data$classification==.env$classification)
  reference_band<-range(wrap12(c(20+27/60,23+55/60)-c(2+24/60,6+52/60)))
  for(axis in c('clock','phase')) {
    z<-filter(zall,.data$axis==.env$axis) |> mutate(context=factor(context,c('Overall','Work','Free')))
    lab<-distinct(z,context,group,participants,participant_days) |>
      mutate(text=paste0(chronotype_plot_label(group),': ',participants,' participants; ',participant_days,' participant-days'),v=if_else(group=='Morning',1.2,2.7))
    p<-ggplot(z,aes(x,estimate_melEDI,color=group,fill=group))
    if(axis=='phase')p<-p+dlmo_reference_band()
    p<-p+geom_ribbon(aes(ymin=lower_melEDI,ymax=upper_melEDI),alpha=.15,color=NA)+geom_line()+facet_wrap(~context,nrow=1)+
      geom_text(data=lab,aes(x=Inf,y=Inf,label=text,color=group,vjust=v),inherit.aes=FALSE,hjust=1.02,size=2.6,show.legend=FALSE)+
      scale_color_manual(values=c(Morning='#0072B2',Evening='#D55E00'))+scale_fill_manual(values=c(Morning='#0072B2',Evening='#D55E00'))+
      time_scale(axis)+intensity_scale(limits=range(c(0,zall$lower_melEDI,zall$upper_melEDI)),expand=expansion(mult=c(.03,.25)))+
      labs(x=axis_title(axis),y='Melanopic EDI (lx; symlog)',color=NULL,fill=NULL)+theme(legend.position='none')
    p<-p+labs(tag=if(axis=='clock')'A'else'B')+theme(plot.tag=element_text(face='bold'))
    panels<-c(panels,list(p))
  }
  patchwork::wrap_plots(panels,ncol=1)+
    patchwork::plot_annotation(title=paste(classification,'sensitivity: near-eye exposure'),
      subtitle='68 participants; 379 participant-days; 17,558 half-hours. Common photoperiod and day-type weights; site/person/day terms excluded.',
      caption='A: clock time. B: circadian phase. Morning/evening pools definite and moderate MEQ categories; MSFsc uses the same people and original median.\nBands: 95% pointwise coefficient intervals. Lux values are back-transformed log-scale means; symlog is linear below 1 lx.\nStriped band: DLMO reference. Time-specific support remains in the exported tables.')
}

meq_factor_figure <- function(contrasts) {
  panels<-list()
  for(axis in c('clock','phase')) {
    z<-filter(contrasts,.data$axis==.env$axis) |> mutate(context=factor(context,c('Overall','Work','Free')))
    lab<-distinct(z,classification,context,evening_participants,evening_days,morning_participants,morning_days) |>
      mutate(text=paste0(classification,': evening ',evening_participants,' / ',evening_days,'; morning ',morning_participants,' / ',morning_days),
        v=if_else(classification=='MEQ',1.2,2.7))
    p<-ggplot(z,aes(x,ratio,color=classification,fill=classification))
    if(axis=='phase')p<-p+dlmo_reference_band()
    p<-p+geom_hline(yintercept=1,color='grey50')+
      geom_ribbon(aes(ymin=ratio_lower,ymax=ratio_upper),alpha=.12,color=NA)+geom_line()+facet_wrap(~context,nrow=1)+
      geom_text(data=lab,aes(x=Inf,y=Inf,label=text,color=classification,vjust=v),inherit.aes=FALSE,hjust=1.02,size=2.8,show.legend=FALSE)+
      scale_color_manual(values=c(MEQ='#009E73',MSFsc='#CC79A7'))+scale_fill_manual(values=c(MEQ='#009E73',MSFsc='#CC79A7'))+
      time_scale(axis)+factor_scale(limits=range(c(1,contrasts$ratio_lower,contrasts$ratio_upper)),expand=expansion(mult=c(.04,.24)))+
      labs(x=axis_title(axis),y='Evening / morning factor (melEDI + 0.1)',color=NULL,fill=NULL,tag=if(axis=='clock')'A'else'B')+
      theme(plot.tag=element_text(face='bold'))
    panels[[axis]]<-p
  }
  patchwork::wrap_plots(panels,ncol=1,guides='collect')+
    patchwork::plot_annotation(title='MEQ and MSFsc grouping on the same restricted sample',
      subtitle='68 participants; 379 participant-days; 17,558 half-hours. Counts: participants / participant-days by chronotype direction.',
      caption='A: clock time. B: circadian phase. Striped band: DLMO reference. MSFsc uses the original median; MEQ pools moderate and definite types.\nFactors refer to melanopic EDI + 0.1. Bands cover all three day contexts and times at 95% within each coordinate and classification.\nDescriptive sensitivity estimates under the Gaussian reference variance model; they do not test equivalence of the classifications.')
}

run_meq <- function(prepared=readRDS('results/data/prepared.rds')) {
  raw<-meq_prepare(prepared);prep<-fit_preprocess(raw,prepared$cutoffs)
  save_table(bind_cols(tibble(MSFsc_median=prepared$cutoffs[2],k_population=12L,k_participant=10L,
    covariance_draws=cfg$draws,seed=cfg$seed,response='log10(melanopic EDI + 0.1)',
    formula=paste(deparse(model_formula()),collapse=' ')),count_sample(raw)),'meq_analysis_contract','results/provenance')
  curves<-contrasts<-windows<-summaries<-verdicts<-list();i<-0L
  for(classification in c('MEQ','MSFsc'))for(axis in c('clock','phase')) {
    i<-i+1L;b<-meq_fit(raw,prep,axis,classification)
    summaries[[i]]<-mutate(model_summary(b),classification=classification)
    di<-diagnose_model(b);verdicts[[i]]<-mutate(di$flags,classification=classification)
    inf<-meq_inference(b);curves[[i]]<-inf$curves;contrasts[[i]]<-inf$contrasts;windows[[i]]<-inf$windows
    rm(b,di,inf);gc()
  }
  curves<-bind_rows(curves);contrasts<-bind_rows(contrasts);windows<-bind_rows(windows)
  save_table(bind_rows(summaries),'meq_model_fit');save_table(bind_rows(verdicts),'meq_diagnostic_verdicts')
  save_table(curves,'meq_all_curves');save_table(contrasts,'meq_all_contrasts');save_table(windows,'meq_candidate_windows')
  for(classification in c('MEQ','MSFsc'))save_plot(meq_pattern_figure(curves,classification),paste0('meq_patterns_',tolower(classification)),14,11)
  save_plot(meq_factor_figure(contrasts),'meq_factor_comparison',14,9)
  meq_verify(raw,curves,contrasts,prepared$cutoffs)
  writeLines(c(capture.output(sessionInfo()),'','Command: Rscript --vanilla R/meq.R',
    paste('OpenMP preload:',omp),
    paste('Input prepared.rds MD5:',tools::md5sum('results/data/prepared.rds'))),'results/provenance/meq_session.txt')
  log_step('Completed MEQ sensitivity')
  invisible(list(curves=curves,contrasts=contrasts,windows=windows))
}

meq_verify <- function(raw,curves,contrasts,cutoffs) {
  # Equality checks use the analytical registry and the actual model frames.
  checks<-list();check<-function(name,ok){stopifnot(ok);checks[[length(checks)+1L]]<<-tibble(check=name,passed=isTRUE(ok))}
  for(classification in c('meq_defined','meq_matched_msfsc')) {
    a<-readRDS(paste0('results/models/',classification,'_clock.rds'));b<-readRDS(paste0('results/models/',classification,'_phase.rds'))
    check(paste(classification,'clock/phase observations and AR order match'),identical(a$data$participant_day_key,b$data$participant_day_key)&&
      identical(a$data$clock_bin,b$data$clock_bin)&&identical(a$data$response,b$data$response)&&identical(a$data$AR_start,b$data$AR_start))
    check(paste(classification,'no link crosses a day or gap'),all(a$data$AR_start[!duplicated(a$data$participant_day_key)])&&all(a$data$elapsed_from_previous_seconds[!a$data$AR_start]==0))
    rm(a,b);gc()
  }
  for(axis in c('clock','phase')) {
    a<-readRDS(paste0('results/models/meq_defined_',axis,'.rds'));b<-readRDS(paste0('results/models/meq_matched_msfsc_',axis,'.rds'))
    check(paste(axis,'paired classifications use identical ordered observations and photoperiod centers'),
      identical(a$data$participant_day_key,b$data$participant_day_key)&&identical(a$data$clock_bin,b$data$clock_bin)&&
      identical(a$data$response,b$data$response)&&identical(a$data$AR_start,b$data$AR_start)&&identical(a$data$photo,b$data$photo))
    rm(a,b);gc()
  }
  check('Factor signs retain coefficient-scale direction',all((contrasts$ratio>1)==(contrasts$estimate>0)))
  check('Candidate windows obey local support rule',!any(contrasts$candidate_interval & contrasts$minimum_cell_participants<cfg$support_n))
  check('Primary MEQ sample',identical(as.integer(unlist(count_sample(raw))),c(68L,379L,17558L)))
  check('One observed half-hour per participant-day and clock bin',!anyDuplicated(raw[c('participant_day_key','clock_bin')]))
  check('Same-subset MSFsc groups retain frozen median',all(raw$msf_direction==ifelse(raw$msf_hours<=cutoffs[2],'Morning','Evening')))
  check('MSFsc-relative display wraps correctly',isTRUE(all.equal(raw$phase_display,wrap12(raw$clock-raw$msf_hours))))
  for(classification in c('MEQ','MSFsc'))for(axis in c('clock','phase'))for(context in c('Overall','Work','Free'))for(group in c('Morning','Evening')) {
    d<-raw;d$quartile<-if(classification=='MEQ')d$meq_direction else d$msf_direction
    z<-filter(curves,.data$classification==.env$classification,.data$axis==.env$axis,.data$context==.env$context,.data$group==.env$group)
    s<-meq_support(d,group,context,axis,z$x)
    check(paste(classification,axis,context,group,'curve totals and local support match registry'),
      all(z$participants==s$total$participants)&all(z$participant_days==s$total$participant_days)&all(z$half_hours==s$total$half_hours)&
        identical(as.integer(z$support_participants),as.integer(s$n))&identical(as.integer(z$support_days),as.integer(s$days)))
    ct<-filter(contrasts,.data$classification==.env$classification,.data$axis==.env$axis,.data$context==.env$context)
    prefix<-tolower(group)
    check(paste(classification,axis,context,group,'contrast counts match registry'),
      all(ct[[paste0(prefix,'_participants')]]==s$total$participants)&all(ct[[paste0(prefix,'_days')]]==s$total$participant_days)&
        all(ct[[paste0(prefix,'_half_hours')]]==s$total$half_hours))
  }
  save_table(bind_rows(checks),'meq_acceptance','results/provenance')
  invisible(bind_rows(checks))
}

if(sys.nframe()==0L) {
  source('R/setup.R')
  for(f in c('data','models','inference','diagnostics','plots'))source(paste0('R/',f,'.R'))
  cfg$nthreads<-1L
  run_meq()
}
