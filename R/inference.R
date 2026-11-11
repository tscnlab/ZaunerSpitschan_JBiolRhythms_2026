cell_registry <- function() expand_grid(quartile=paste0('Q',1:4),daytype=c('Work','Free')) |>
  mutate(cell=paste(quartile,daytype,sep='_'),direction=if_else(quartile %in% c('Q1','Q2'),'Morning','Evening'),
    extremity=if_else(quartile %in% c('Q2','Q3'),'Moderate','Extreme'))
contrast_registry <- function(free_weight) {
  cells <- cell_registry(); out <- list(); weights <- list()
  add <- function(id,family,context,description,values) {
    out[[id]] <<- tibble(contrast=id,family=family,context=context,description=description)
    weights[[id]] <<- values
  }
  q_weight <- function(q,context) {
    w <- as.numeric(cells$quartile==q)
    if(context=='Overall') w <- w*ifelse(cells$daytype=='Free',free_weight,1-free_weight)
    else w <- w*as.numeric(cells$daytype==context)
    w
  }
  for(context in c('Overall','Work','Free')) {
    for(severity in c('Moderate','Extreme')) {
      q <- if(severity=='Moderate')c('Q2','Q3')else c('Q1','Q4')
      add(paste('evening_morning',severity,context,sep='_'),'Evening minus morning',context,paste('Evening minus morning:',severity),q_weight(q[2],context)-q_weight(q[1],context))
    }
    for(direction in c('Morning','Evening')) {
      q <- if(direction=='Morning')c('Q2','Q1')else c('Q3','Q4')
      add(paste('extreme_moderate',direction,context,sep='_'),'Extreme minus moderate',context,paste('Extreme minus moderate:',direction),q_weight(q[2],context)-q_weight(q[1],context))
    }
  }
  for(q in paste0('Q',1:4)) add(paste0('free_work_',q),'Free minus work','Free versus work',paste('Free minus work:',q),q_weight(q,'Free')-q_weight(q,'Work'))
  de <- function(context) q_weight('Q4',context)-q_weight('Q1',context)-q_weight('Q3',context)+q_weight('Q2',context)
  dd <- function(context) (q_weight('Q4',context)+q_weight('Q3',context)-q_weight('Q1',context)-q_weight('Q2',context))/2
  ed <- function(context) (q_weight('Q4',context)+q_weight('Q1',context)-q_weight('Q3',context)-q_weight('Q2',context))/2
  add('direction_extremity','Interactions','Overall','Direction by extremity',de('Overall'))
  add('direction_daytype','Interactions','Free versus work','Direction by day type',dd('Free')-dd('Work'))
  add('extremity_daytype','Interactions','Free versus work','Extremity by day type',ed('Free')-ed('Work'))
  add('direction_extremity_daytype','Interactions','Free versus work','Direction by extremity by day type',de('Free')-de('Work'))
  W <- do.call(rbind,weights); colnames(W)<-cells$cell
  stopifnot(all(abs(rowSums(W))<1e-12),nrow(W)==20)
  list(table=bind_rows(out),W=W,cells=cells)
}
standardized_cells <- function(bundle,reference=NULL) {
  ref <- if(is.null(reference))reference_days(bundle$data)else reference
  ref <- ref |> summarise(w=sum(w),.by=photo)
  x <- if(bundle$axis=='clock')seq(0,24-cfg$step,by=cfg$step)else seq(-12,12-cfg$step,by=cfg$step)
  t <- wrap24(x); n <- length(t)
  nd <- prediction_template(bundle,n);nd$t<-t;nd$photo<-0
  base <- population_matrix(bundle,nd)
  # Exact weighted integration of photoperiod basis columns. Other columns do not depend on photo.
  for(sm in bundle$fit$smooth) if('photo' %in% sm$term) {
    ix <- sm$first.para:sm$last.para; avg <- matrix(0,n,length(ix))
    for(j in seq_len(nrow(ref))) { nd$photo<-ref$photo[j]; avg<-avg + ref$w[j]*mgcv::PredictMat(sm,nd) }
    base[,ix] <- avg
  }
  nd$photo<-0
  Xzero<-population_matrix(bundle,nd)
  delta<-base-Xzero
  cells <- cell_registry()
  matrices <- lapply(seq_len(nrow(cells)),function(j){
    z<-nd;for(v in c('quartile','direction','extremity','daytype'))z[[v]]<-cells[[v]][j]
    population_matrix(bundle,z)+delta
  })
  names(matrices)<-cells$cell
  list(x=x,X=matrices,cells=cells,reference=ref)
}
counts_for_cells <- function(d,cells) {
  sub<-filter(d,paste(quartile,daytype,sep='_') %in% cells)
  count_sample(sub)
}
support_for_cells <- function(d,cells,axis,x) {
  sub<-filter(d,paste(quartile,daytype,sep='_') %in% cells)
  origin<-if(axis=='clock')0 else -12
  coord<-if(axis=='clock')sub$clock else sub$phase_display
  sub$bin<-floor((coord-origin)/.5)
  local<-sub |> summarise(n=n_distinct(participant_key),days=n_distinct(participant_day_key),.by=bin)
  i<-match(floor((x-origin)/.5),local$bin)
  tibble(n=replace_na(local$n[i],0L),days=replace_na(local$days[i],0L))
}
covariance_info <- function(f) {
  unconditional<-!is.null(f$Vc)
  list(V=vcov(f,unconditional=unconditional),unconditional=unconditional,
       label=if(unconditional)'Bayesian coefficient covariance with smoothing-parameter correction (Vc)'else'Bayesian coefficient covariance conditional on smoothing parameters (Vp)')
}
linear_summary <- function(X,beta,V,critical=1.95996398454) {
  keep<-which(colSums(abs(X))>1e-12); xx<-X[,keep,drop=FALSE]
  est<-drop(xx%*%beta[keep]);se<-sqrt(pmax(0,rowSums((xx%*%V[keep,keep,drop=FALSE])*xx)))
  tibble(estimate=est,se=se,lower=est-critical*se,upper=est+critical*se)
}
infer_curves <- function(bundle,reference=NULL,axis_coverage=cfg$curve_coverage,export=TRUE,
  pattern_coverage=if(bundle$label=='primary')cfg$curve_coverage else NULL) {
  stopifnot(length(axis_coverage)==1L,is.finite(axis_coverage),axis_coverage>.5,axis_coverage<1)
  log_step('Inference',bundle$label,bundle$axis)
  cells<-standardized_cells(bundle,reference);d<-bundle$data; f<-bundle$fit
  day_ref<-reference_days(d);free_weight<-sum(day_ref$w[day_ref$daytype=='Free'])
  cr<-contrast_registry(free_weight);V<-covariance_info(f);beta<-coef(f)
  curves<-list();pattern_matrices<-list()
  for(context in c('Overall','Work','Free')) for(q in paste0('Q',1:4)) {
    cc<-if(context=='Overall')paste(q,c('Work','Free'),sep='_')else paste(q,context,sep='_')
    X<-if(context=='Overall')(1-free_weight)*cells$X[[cc[1]]]+free_weight*cells$X[[cc[2]]]else cells$X[[cc]]
    nr<-counts_for_cells(d,cc);su<-support_for_cells(d,cc,bundle$axis,cells$x)
    key<-paste(q,context,sep='_');pattern_matrices[[key]]<-X
    curves[[key]]<-bind_cols(tibble(x=cells$x,quartile=q,context=context,
      direction=if(q %in% c('Q1','Q2'))'Morning'else'Evening',extremity=if(q %in% c('Q2','Q3'))'Moderate'else'Extreme'),
      linear_summary(X,beta,V$V)) |> mutate(participants=nr$participants,participant_days=nr$participant_days,half_hours=nr$half_hours,support_n=su$n,support_days=su$days)
  }
  D<-lapply(seq_len(nrow(cr$W)),function(i) Reduce(`+`,Map(function(X,w)X*w,cells$X,cr$W[i,])))
  names(D)<-rownames(cr$W)
  DD<-do.call(rbind,D);keep<-which(colSums(abs(DD))>1e-12)
  ds<-linear_summary(DD,beta,V$V)
  set.seed(cfg$seed+if(bundle$axis=='phase')1L else 0L)
  # gratia uses the installed mgcv covariance, including Vc when available.
  draws<-as.matrix(gratia::gaussian_draws(f,n=cfg$draws,unconditional=V$unconditional))
  centered<-sweep(draws[,keep,drop=FALSE],2,beta[keep],'-')
  maxz<-numeric(nrow(draws));usable<-ds$se>1e-10
  for(ii in split(seq_len(nrow(draws)),ceiling(seq_len(nrow(draws))/250))) {
    zz<-centered[ii,,drop=FALSE]%*%t(DD[usable,keep,drop=FALSE])
    zz<-sweep(zz,2,ds$se[usable],'/');maxz[ii]<-apply(abs(zz),1,max)
  }
  # Each time coordinate addresses a separate question. Cover the full
  # registered contrast family and grid within that coordinate at 95%.
  critical<-unname(quantile(maxz,axis_coverage,type=8))
  # Primary exposure bands cover all four quartiles and three day contexts
  # across the grid within this coordinate. Use the same joint coefficient
  # draws as the contrast family, without an adjustment across coordinates.
  pattern_critical<-qnorm(.975);pattern_maxz<-NULL
  if(!is.null(pattern_coverage)) {
    stopifnot(length(pattern_coverage)==1L,pattern_coverage>.5,pattern_coverage<1)
    XX<-do.call(rbind,pattern_matrices);pp<-which(colSums(abs(XX))>1e-12)
    ps<-linear_summary(XX,beta,V$V);usable_p<-ps$se>1e-10
    pattern_maxz<-numeric(nrow(draws))
    for(ii in split(seq_len(nrow(draws)),ceiling(seq_len(nrow(draws))/250))) {
      zz<-sweep(draws[ii,pp,drop=FALSE],2,beta[pp],'-')%*%t(XX[usable_p,pp,drop=FALSE])
      zz<-sweep(zz,2,ps$se[usable_p],'/');pattern_maxz[ii]<-apply(abs(zz),1,max)
    }
    pattern_critical<-unname(quantile(pattern_maxz,pattern_coverage,type=8))
  }
  curves<-bind_rows(curves) |> mutate(lower=estimate-pattern_critical*se,
    upper=estimate+pattern_critical*se,
    band_type=if(is.null(pattern_coverage))'pointwise'else'simultaneous',
    band_coverage=if(is.null(pattern_coverage)).95 else pattern_coverage,
    band_critical=pattern_critical)
  contrasts<-list()
  for(i in seq_len(nrow(cr$W))) {
    w<-cr$W[i,];pos<-names(w)[w>0];neg<-names(w)[w<0]
    a<-counts_for_cells(d,pos);b<-counts_for_cells(d,neg)
    support<-vapply(names(w)[w!=0],function(cc)support_for_cells(d,cc,bundle$axis,cells$x)$n,numeric(length(cells$x)))
    supported<-apply(support,1,min)>=cfg$support_n
    z<-linear_summary(D[[i]],beta,V$V,critical)
    contrasts[[i]]<-bind_cols(cr$table[i,][rep(1,length(cells$x)),],tibble(x=cells$x),z) |>
      mutate(supported=supported,minimum_cell_participants=apply(support,1,min),
        positive_participants=a$participants,positive_days=a$participant_days,positive_half_hours=a$half_hours,
        negative_participants=b$participants,negative_days=b$participant_days,negative_half_hours=b$half_hours,
        ratio=10^estimate,ratio_lower=10^lower,ratio_upper=10^upper,
        significant=supported & (lower>0|upper<0))
  }
  result<-list(curves=curves,contrasts=bind_rows(contrasts),registry=cr,cell_design=cells,
    pattern_matrices=pattern_matrices,critical=critical,covariance=V$label,free_weight=free_weight,
    individual_draws=draws[seq_len(1000),,drop=FALSE],axis_coverage=axis_coverage,
    minimum_support=cfg$support_n,simulated_maxima=maxz,
    pattern_coverage=pattern_coverage,pattern_critical=pattern_critical,
    pattern_simulated_maxima=pattern_maxz)
  result$windows<-significant_windows(result,bundle)
  if(!export)return(result)
  stem<-paste(bundle$label,bundle$axis,sep='_')
  save_table(intensity_display(result$curves),paste0(stem,'_curves'));save_table(contrast_display(result$contrasts),paste0(stem,'_contrasts'))
  save_table(window_display(result$windows),paste0(stem,'_significant_windows'))
  save_table(cbind(cr$table,as.data.frame(cr$W)),paste0(stem,'_contrast_registry'))
  save_table(bind_cols(tibble(model=bundle$label,axis=bundle$axis,covariance=V$label,draws=cfg$draws,critical=critical,
    axis_family_coverage=axis_coverage,
    inference_family='All 20 registered contrasts and grid points within this time coordinate',
    across_coordinate_adjustment='None; clock time and phase address separate questions',
    exposure_band_type=if(is.null(pattern_coverage))'pointwise'else'simultaneous',
    exposure_band_coverage=if(is.null(pattern_coverage)).95 else pattern_coverage,
    exposure_band_critical=pattern_critical,
    exposure_band_family=if(is.null(pattern_coverage))'Each grid point separately'else'All 12 quartile-by-context curves and grid points within this time coordinate',
    minimum_cell_support=cfg$support_n,free_day_weight=free_weight),count_sample(d)),paste0(stem,'_inference_metadata'))
  result
}

# Average the design vectors before calculating uncertainty, retaining covariance
# across phase points. The selected phase windows are treated as fixed here.
interval_average_contrasts <- function(inference,bundle,
  contrast='evening_morning_Extreme_Overall',coverage=.95) {
  stopifnot(contrast %in% rownames(inference$registry$W),coverage>0,coverage<1)
  w<-inference$registry$W[contrast,]
  D<-Reduce(`+`,Map(function(X,w)X*w,inference$cell_design$X,w))
  windows<-filter(inference$windows,.data$contrast==.env$contrast) |> arrange(start)
  V<-covariance_info(bundle$fit)
  out<-lapply(seq_len(nrow(windows)),function(i) {
    ix<-which(wrap24(inference$cell_design$x-windows$start[i])<windows$duration_hours[i]-1e-10)
    stopifnot(length(ix)>0)
    X<-matrix(colMeans(D[ix,,drop=FALSE]),nrow=1)
    z<-linear_summary(X,coef(bundle$fit),V$V,qnorm((1+coverage)/2))
    bind_cols(windows[i,c('axis','contrast','context','start','end','duration_hours')],z) |>
      mutate(ratio=10^estimate,ratio_lower=10^lower,ratio_upper=10^upper,
        coverage=coverage,grid_points=length(ix),covariance=V$label,
        interval_method='Mean log contrast; selected window fixed; no selection or multiplicity adjustment')
  })
  bind_rows(out)
}
significant_windows <- function(inference,bundle) {
  d<-bundle$data;axis<-bundle$axis;origin<-if(axis=='clock')0 else -12
  out<-list()
  for(id in inference$registry$table$contrast) {
    z<-filter(inference$contrasts,contrast==id) |> arrange(x)
    sign<-ifelse(z$significant,ifelse(z$estimate>0,1L,-1L),0L)
    rr<-rle(sign);ends<-cumsum(rr$lengths);starts<-ends-rr$lengths+1L
    groups<-Map(seq.int,starts,ends);values<-rr$values
    if(length(groups)>1L && values[1]!=0 && values[1]==tail(values,1)) {
      groups[[1]]<-c(tail(groups,1)[[1]],groups[[1]]);groups<-head(groups,-1);values<-head(values,-1)
    }
    w<-inference$registry$W[id,]
    for(j in which(values!=0)) {
      ix<-groups[[j]];start<-z$x[ix[1]]-cfg$step/2;duration<-length(ix)*cfg$step
      start<-wrap24(start-origin)+origin;end<-wrap24(start+duration-origin)+origin
      coord<-if(axis=='clock')d$clock else d$phase_display
      inside<-wrap24(coord-start)<duration-1e-10
      a<-counts_for_cells(d[inside,],names(w)[w>0]);b<-counts_for_cells(d[inside,],names(w)[w<0])
      if(min(a$participants,b$participants)<cfg$support_n)next
      out[[length(out)+1L]]<-tibble(model=bundle$label,axis=axis,contrast=id,description=z$description[1],context=z$context[1],
        direction=if(values[j]>0)'Positive'else'Negative',start=start,end=end,wraps_boundary=start>=end,duration_hours=duration,
        mean_difference=mean(z$estimate[ix]),minimum_difference=min(z$estimate[ix]),maximum_difference=max(z$estimate[ix]),
        mean_ratio=10^mean(z$estimate[ix]),positive_participants=a$participants,positive_days=a$participant_days,positive_half_hours=a$half_hours,
        negative_participants=b$participants,negative_days=b$participant_days,negative_half_hours=b$half_hours)
    }
  }
  if(length(out))bind_rows(out)else tibble(model=character(),axis=character(),contrast=character(),description=character(),context=character(),direction=character(),start=double(),end=double(),wraps_boundary=logical(),duration_hours=double(),mean_difference=double(),minimum_difference=double(),maximum_difference=double(),mean_ratio=double(),positive_participants=integer(),positive_days=integer(),positive_half_hours=integer(),negative_participants=integer(),negative_days=integer(),negative_half_hours=integer())
}
