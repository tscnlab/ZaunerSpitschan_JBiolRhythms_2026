blocks <- c('direction','extremity','daytype','site','photoperiod','participant','participant_day')
model_formula <- function(included = blocks, k = cfg$k, kp = cfg$kp) {
  # mgcv 1.9-4 misstates the rank for linked cyclic sz penalties (id != NULL).
  # Its standard separate-penalty constructor fits the requested cyclic basis.
  sz <- function(f) sprintf("s(t,%s,bs='sz',xt=list(bs='cc'),k=%d)",paste(f,collapse=','),k)
  terms <- sprintf("s(t,bs='cc',k=%d)",k)
  grouping <- if (all(c('direction','extremity') %in% included)) 'quartile' else intersect(c('direction','extremity'),included)
  if(length(grouping)) terms <- c(terms,sz(grouping))
  if('daytype' %in% included) {
    terms <- c(terms,sz('daytype'))
    if(length(grouping)) terms <- c(terms,sz(c(grouping,'daytype')))
  }
  if('site' %in% included) terms <- c(terms,sz('site'))
  if('photoperiod' %in% included) terms <- c(terms,"s(photo,bs='cr',k=5)",sprintf("ti(t,photo,bs=c('cc','cr'),k=c(%d,5))",k))
  if('participant' %in% included) terms <- c(terms,sprintf("s(t,participant,bs='fs',xt=list(bs='cc'),k=%d)",kp))
  if('participant_day' %in% included) terms <- c(terms,"s(participant_day,bs='re')")
  as.formula(paste('response ~',paste(terms,collapse=' + ')))
}
lag_cor <- function(e, start, lag = 1L) {
  i <- seq.int(lag+1L,length(e)); j <- i-lag; s <- cumsum(start); ok <- s[i]==s[j] & is.finite(e[i]) & is.finite(e[j])
  if(sum(ok)<3) return(NA_real_)
  cor(e[i[ok]],e[j[ok]])
}
fit_bam <- function(d, formula = model_formula(), rho = 0, label = '', weights = NULL) {
  log_step('Fit',label,'n =',nrow(d),'rho =',round(rho,4))
  warnings <- character()
  d$.analysis_weight <- if(is.null(weights))rep(1,nrow(d))else weights
  tm <- system.time(f <- withCallingHandlers(bam(formula, data=d, family=gaussian(), method='fREML', discrete=TRUE,
    knots=list(t=c(0,24)),rho=rho,AR.start=AR_start,weights=.analysis_weight,nthreads=cfg$nthreads,gc.level=1,drop.unused.levels=TRUE),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}))
  attr(f,'analysis_warnings') <- unique(warnings); attr(f,'elapsed_seconds') <- unname(tm[['elapsed']])
  log_step('Finished',label,'in',round(tm[['elapsed']],1),'seconds')
  f
}
fit_pair <- function(raw, prep, axis, label, k=cfg$k, kp=cfg$kp) {
  d <- apply_preprocess(raw,prep,axis); formula <- model_formula(k=k,kp=kp)
  signature <- digest::digest(list(d,formula=deparse(formula),cfg$version,mgcv=as.character(packageVersion('mgcv'))))
  path <- file.path('results/models',paste0(label,'_',axis,'.rds'))
  if(file.exists(path)) { old <- readRDS(path); if(identical(old$signature,signature)) return(old) }
  preliminary <- fit_bam(d,formula,0,paste(label,axis,'preliminary'))
  set.seed(cfg$seed)
  preliminary_check <- mgcv::k.check(preliminary,subsample=5000,n.rep=100)
  preliminary_patterns <- tibble(t=d$t,residual=residuals(preliminary),fitted=fitted(preliminary),quartile=d$quartile,daytype=d$daytype) |>
    mutate(hour=floor(t)) |> summarise(mean_residual=mean(residual),sd_residual=sd(residual),half_hours=n(),.by=c(hour,quartile,daytype))
  # Preserve the preliminary mean-structure checks before estimating dependence.
  rho <- pmin(.95,pmax(-.95,lag_cor(residuals(preliminary),d$AR_start)))
  save_table(cbind(term=rownames(preliminary_check),as.data.frame(preliminary_check)),paste0(label,'_',axis,'_preliminary_k'),'results/diagnostics')
  save_table(preliminary_patterns,paste0(label,'_',axis,'_preliminary_patterns'),'results/diagnostics')
  rm(preliminary); gc()
  f <- fit_bam(d,formula,rho,paste(label,axis,'AR1'))
  out <- list(fit=f,data=d,prep=prep,axis=axis,label=label,rho=rho,signature=signature,preliminary_k=preliminary_check)
  atomic_rds(out,path)
  out
}
model_summary <- function(bundle) {
  # These summaries do not use random-effect significance tests.
  f<-bundle$fit; s<-summary(f,re.test=FALSE)
  bind_cols(tibble(model=bundle$label,axis=bundle$axis,rho=bundle$rho,
    adjusted_R2=s$r.sq,deviance_explained=s$dev.expl,
    ordinary_R2=1-sum(residuals(f)^2)/sum((bundle$data$response-mean(bundle$data$response))^2),
    edf=sum(f$edf),coefficients=length(coef(f)),rank=f$rank,
    convergence=if(!is.null(f$outer.info$conv))as.character(f$outer.info$conv) else as.character(f$converged),
    warnings=paste(attr(f,'analysis_warnings'),collapse='; ')),count_sample(bundle$data))
}
exclude_labels <- function(f, site=TRUE, participant=TRUE, day=TRUE) {
  labs <- vapply(f$smooth,function(s)s$label,character(1))
  labs[(site & grepl(',site)',labs,fixed=TRUE)) | (participant & grepl(',participant)',labs,fixed=TRUE)) | (day & labs=='s(participant_day)')]
}
prediction_template <- function(bundle, n=1L) bundle$data[rep(1L,n),,drop=FALSE]
set_factors <- function(nd,bundle) {
  for(v in c('quartile','direction','extremity','daytype','site','participant','participant_day')) nd[[v]] <- factor(nd[[v]],levels=levels(bundle$data[[v]]))
  nd
}
population_matrix <- function(bundle, nd) {
  predict(bundle$fit,newdata=set_factors(nd,bundle),type='lpmatrix',exclude=exclude_labels(bundle$fit),discrete=FALSE)
}
reference_days <- function(d) d |> distinct(participant_key,participant_day_key,photo,daytype) |>
  group_by(participant_key) |> mutate(w=1/n()) |> ungroup() |> mutate(w=w/n_distinct(participant_key))
