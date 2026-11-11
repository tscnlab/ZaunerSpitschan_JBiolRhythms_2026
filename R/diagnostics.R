qr_concurvity <- function(f) {
  # Same unweighted design and projection definition, with a rank-tolerant QR.
  # Avoids the triangular inverse/SVD that fails for some aliased fs/sz designs.
  X<-qr.R(qr(model.matrix(f),tol=0,LAPACK=FALSE))
  indices<-lapply(f$smooth,function(s)s$first.para:s$last.para)
  names(indices)<-vapply(f$smooth,function(s)s$label,character(1))
  indices<-c(list(para=seq_len(min(unlist(indices))-1L)),indices)
  bind_rows(lapply(names(indices),function(term){
    j<-indices[[term]];Xi<-X[,-j,drop=FALSE];Xj<-X[,j,drop=FALSE]
    q<-qr(Xi,tol=1e-8,LAPACK=FALSE);Q<-qr.Q(q)[,seq_len(q$rank),drop=FALSE]
    projected<-Q%*%crossprod(Q,Xj);beta<-coef(f)[j]
    tibble(.type=c('observed','estimate'),.term=term,
      .concurvity=c(sum((projected%*%beta)^2)/sum((Xj%*%beta)^2),sum(projected^2)/sum(Xj^2)))
  }))
}
diagnose_model <- function(bundle) {
  f<-bundle$fit;d<-bundle$data;stem<-paste(bundle$label,bundle$axis,sep='_')
  log_step('Diagnostics',stem)
  set.seed(cfg$seed+201L)
  kc<-mgcv::k.check(f,subsample=5000,n.rep=200)
  kt<-cbind(term=rownames(kc),as.data.frame(kc));rownames(kt)<-NULL
  kt<-bind_cols(as_tibble(kt),count_sample(d)[rep(1,nrow(kt)),])
  save_table(kt,paste0(stem,'_basis'),'results/diagnostics')
  capture.output({pdf(paste0('results/figures/',stem,'_gam_check.pdf'),width=9,height=8);mgcv::gam.check(f,k.rep=100);dev.off()},file=paste0('results/diagnostics/',stem,'_gam_check.txt'))
  ap<-try(gratia::appraise(f),silent=TRUE)
  if(!inherits(ap,'try-error'))save_plot(ap+patchwork::plot_annotation(title=paste(stem,'model checks'),subtitle=sample_label(d)),paste0(stem,'_appraise'),11,8)
  cv<-try(gratia::model_concurvity(f),silent=TRUE)
  if(inherits(cv,'try-error')) {
    writeLines(c(as.character(cv),'Fallback: unweighted model-design projections using QR with rank tolerance 1e-8. Worst-case SVD index omitted.'),paste0('results/diagnostics/',stem,'_concurvity_fallback.txt'))
    cv<-qr_concurvity(f)
  }
  save_table(bind_cols(as.data.frame(cv),count_sample(d)[rep(1,nrow(cv)),]),paste0(stem,'_concurvity'),'results/diagnostics')
  e<-residuals(f);std<-f$std.rsd
  stopifnot(length(e)==nrow(d),length(std)==nrow(d))
  acf_one<-function(z,lag){if(length(z)<=lag+2||sd(z)==0)return(NA_real_);cor(z[(lag+1):length(z)],z[1:(length(z)-lag)])}
  series<-split(seq_len(nrow(d)),d$true_elapsed_sequence_id)
  acf<-bind_rows(lapply(names(series),function(id){ix<-series[[id]];tibble(sequence=id,participant_key=d$participant_key[ix[1]],participant_day_key=d$participant_day_key[ix[1]],
    lag=1:12,raw=vapply(1:12,function(k)acf_one(e[ix],k),numeric(1)),standardized=vapply(1:12,function(k)acf_one(std[ix],k),numeric(1)),half_hours=length(ix))}))
  save_table(acf,paste0(stem,'_series_acf'),'results/diagnostics')
  acfs<-acf |> pivot_longer(c(raw,standardized),names_to='residual_type',values_to='correlation') |>
    filter(is.finite(correlation)) |> summarise(median=median(correlation),lower=quantile(correlation,.1),upper=quantile(correlation,.9),
      sequences=n_distinct(sequence),participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),.by=c(lag,residual_type))
  save_table(acfs,paste0(stem,'_acf_summary'),'results/diagnostics')
  p<-ggplot(acfs,aes(lag*.5,median,color=residual_type,fill=residual_type))+geom_hline(yintercept=0,color='grey70')+
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.15,color=NA)+geom_line()+geom_point()+
    labs(title=paste(stem,'residual autocorrelation'),subtitle=sample_label(d),x='Lag (hours)',y='Within-sequence correlation',
      caption='Median and 10th to 90th percentiles across eligible true-elapsed sequences. No sequence crosses a day, gap or ambiguous daylight-saving link.')
  save_plot(p,paste0(stem,'_acf'),10,5)
  rd<-d |> mutate(residual=e,standardized=std,fitted=fitted(f),hour=floor(t),decile=ntile(fitted,10))
  patterns<-rd |> summarise(mean_residual=mean(residual),sd_residual=sd(residual),mean_standardized=mean(standardized),
    participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(hour,quartile,daytype))
  save_table(patterns,paste0(stem,'_residual_patterns'),'results/diagnostics')
  p<-ggplot(patterns,aes(hour+.5,mean_residual,color=quartile))+geom_hline(yintercept=0,color='grey50')+geom_line()+
    facet_wrap(~daytype)+labs(title=paste(stem,'mean residuals'),subtitle=sample_label(d),x=if(bundle$axis=='clock')'Clock hour'else'Hours after MSFsc (0 to 24)',y='Mean log-scale residual')
  save_plot(p,paste0(stem,'_residual_patterns'),11,5)
  spread<-rd |> summarise(fitted=mean(fitted),residual_sd=sd(residual),standardized_sd=sd(standardized),
    participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=decile)
  save_table(spread,paste0(stem,'_spread'),'results/diagnostics')
  save_plot(ggplot(spread,aes(fitted,residual_sd))+geom_line()+geom_point()+labs(title=paste(stem,'residual spread'),subtitle=sample_label(d),x='Mean fitted log10(melEDI + 0.1)',y='Residual SD'),paste0(stem,'_spread'),8,5)
  # Verify endpoint and first-derivative continuity for all observed participants and group/site combinations.
  nd<-d |> distinct(participant,quartile,daytype,site,.keep_all=TRUE)
  a<-nd;b<-nd;a$t<-0;b$t<-24
  p0<-predict(f,a,discrete=FALSE);p24<-predict(f,b,discrete=FALSE)
  a$t<-1e-4;b$t<-24-1e-4
  der0<-(predict(f,a,discrete=FALSE)-p0)/1e-4;der24<-(p24-predict(f,b,discrete=FALSE))/1e-4
  joins<-bind_cols(tibble(maximum_endpoint_difference=max(abs(p0-p24)),maximum_derivative_difference=max(abs(der0-der24)),evaluated_rows=nrow(nd)),count_sample(d))
  save_table(joins,paste0(stem,'_cyclic_continuity'),'results/diagnostics')
  inf<-rd |> summarise(mean_residual=mean(residual),RMSE=sqrt(mean(residual^2)),maximum_abs_residual=max(abs(residual)),
    participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(site,participant_key)) |> mutate(participants=1L)
  save_table(inf,paste0(stem,'_participant_residuals'),'results/diagnostics')
  standardized_lag1<-lag_cor(std,d$AR_start)
  centered<-e-mean(e);skew<-mean(centered^3)/mean(centered^2)^1.5;kurt<-mean(centered^4)/mean(centered^2)^2-3
  flags<-bind_cols(tibble(model=bundle$label,axis=bundle$axis,
    convergence_ok=isTRUE(f$converged)&&!any(grepl('not converg|step fail|rank deficient',attr(f,'analysis_warnings'),ignore.case=TRUE)),
    basis_flag=any(kc[,'p-value']<.01 & kc[,'k-index']<.9 & kc[,'edf']/kc[,"k'"]>.8,na.rm=TRUE),
    raw_lag1=lag_cor(e,d$AR_start),standardized_lag1=standardized_lag1,
    residual_skewness=skew,residual_excess_kurtosis=kurt,spread_ratio=max(spread$residual_sd)/min(spread$residual_sd),
    maximum_group_hour_mean_residual=max(abs(patterns$mean_residual)),
    cyclic_ok=joins$maximum_endpoint_difference<1e-6 && joins$maximum_derivative_difference<.01,
    residual_dependence_flag=abs(standardized_lag1)>.15,
    distribution_flag=abs(skew)>1||kurt>3,
    spread_flag=max(spread$residual_sd)/min(spread$residual_sd)>2),count_sample(d))
  save_table(flags,paste0(stem,'_verdict'),'results/diagnostics')
  list(flags=flags,basis=kt,acf=acfs,concurvity=cv,spread=spread,continuity=joins)
}
individual_curves <- function(bundle,inference) {
  f<-bundle$fit;d<-bundle$data;V<-covariance_info(f);beta<-coef(f);x<-inference$cell_design$x;tt<-wrap24(x)
  people<-distinct(d,participant_key,participant,quartile,site)
  result<-timings<-list()
  for(j in seq_len(nrow(people))) {
    pp<-people[j,];ix<-d$participant_key==pp$participant_key;dp<-d[ix,]
    q<-as.character(pp$quartile);base<-inference$pattern_matrices[[paste(q,'Overall',sep='_')]]
    nd<-prediction_template(bundle,length(x));nd$t<-tt;nd$participant<-as.character(pp$participant);nd$quartile<-q
    # Add only the participant's deviation to the same standardized group curve.
    full<-predict(f,set_factors(nd,bundle),type='lpmatrix',exclude=exclude_labels(f,site=TRUE,participant=FALSE,day=TRUE),discrete=FALSE)
    pop<-population_matrix(bundle,nd);X<-base+(full-pop)
    z<-linear_summary(X,beta,V$V);group<-drop(base%*%beta);peak<-x[which.max(z$estimate)];group_peak<-x[which.max(group)]
    keep<-which(colSums(abs(X))>1e-12)
    yy<-inference$individual_draws[,keep,drop=FALSE]%*%t(X[,keep,drop=FALSE]);peaks<-x[max.col(yy,ties.method='first')]
    # Joint draws preserve covariance between the individual and its group reference.
    kb<-which(colSums(abs(base))>1e-12);gg<-inference$individual_draws[,kb,drop=FALSE]%*%t(base[,kb,drop=FALSE])
    gp<-x[max.col(gg,ties.method='first')];shift_draw<-wrap12(peaks-gp)
    offsets<-wrap12(peaks-peak);ci<-quantile(offsets,c(.025,.975),type=8);concentration<-Mod(mean(exp(2i*pi*peaks/24)))
    shift<-wrap12(peak-group_peak);si<-quantile(wrap12(shift_draw-shift),c(.025,.975),type=8)
    nr<-count_sample(dp)
    result[[j]]<-bind_cols(tibble(participant_key=pp$participant_key,quartile=q,x=x,group_estimate=group),z) |>
      mutate(participants=1L,participant_days=nr$participant_days,half_hours=nr$half_hours)
    timings[[j]]<-tibble(participant_key=pp$participant_key,quartile=q,peak=peak,group_peak=group_peak,
      peak_lower_offset=ci[1],peak_upper_offset=ci[2],peak_interval_width=ci[2]-ci[1],peak_concentration=concentration,
      shift_hours=shift,shift_lower=shift+si[1],shift_upper=shift+si[2],
      localized=(ci[2]-ci[1])<=4 && concentration>=.8,
      participants=1L,participant_days=nr$participant_days,half_hours=nr$half_hours)
  }
  stem<-paste(bundle$label,bundle$axis,sep='_');curves<-bind_rows(result);timings<-bind_rows(timings)
  save_table(intensity_display(curves),paste0(stem,'_individual_curves'));save_table(timings,paste0(stem,'_individual_timing'))
  # Every participant appears in the paginated PDF, with its own day count.
  individual_curve_pdf(curves,stem,bundle$axis)
  p<-ggplot(timings,aes(reorder(participant_key,shift_hours),shift_hours,color=localized))+geom_hline(yintercept=0,color='grey70')+
    geom_linerange(aes(ymin=shift_lower,ymax=shift_upper))+geom_point()+coord_flip()+facet_wrap(~quartile,scales='free_y',ncol=2)+
    labs(title=paste(stem,'individual peak shifts'),subtitle=sample_label(d),x=NULL,y='Peak difference from group (hours)',
      caption='Circular differences; 95% joint coefficient-draw intervals. Localized requires a peak interval no wider than 4 h and circular concentration at least 0.8.')
  save_plot(p,paste0(stem,'_timing'),11,14)
  list(curves=curves,timings=timings)
}
