subset_blocks <- function(mask) blocks[as.logical(intToBits(as.integer(mask))[seq_along(blocks)])]
shapley_map <- function() {
  M<-matrix(0,128,7,dimnames=list(as.character(0:127),blocks))
  for(j in 1:7)for(mask in 0:127)if(bitwAnd(mask,bitwShiftL(1L,j-1L))==0L) {
    k<-length(subset_blocks(mask));w<-factorial(k)*factorial(6-k)/factorial(7)
    M[as.character(mask),j]<-M[as.character(mask),j]-w
    M[as.character(bitwOr(mask,bitwShiftL(1L,j-1L))),j]<-M[as.character(bitwOr(mask,bitwShiftL(1L,j-1L))),j]+w
  }
  M
}
subset_fit <- function(mask,bundle) {
  path<-file.path('results/subsets',sprintf('%s_%03d.rds',bundle$axis,mask))
  formula<-model_formula(subset_blocks(mask))
  signature<-digest::digest(list(bundle$signature,bundle$rho,deparse(formula)))
  if(file.exists(path)) {old<-readRDS(path);if(identical(old$signature,signature))return(old)}
  f<-if(mask==127L)bundle$fit else fit_bam(bundle$data,formula,bundle$rho,paste('Shapley',bundle$axis,mask,'of 127'))
  y<-bundle$data$response;pr<-fitted(f)
  out<-list(mask=mask,included=subset_blocks(mask),formula=paste(deparse(formula),collapse=' '),signature=signature,
    rho=bundle$rho,coefficients=coef(f),smoothing_parameters=f$sp,edf=f$edf,
    prediction=as.numeric(pr),R2=1-sum((y-pr)^2)/sum((y-mean(y))^2),
    warnings=attr(f,'analysis_warnings'),converged=f$converged,
    outer_info=f$outer.info,method=f$method,elapsed=attr(f,'elapsed_seconds'))
  atomic_rds(out,path);out
}
hierarchical_day_weights <- function(days,B) {
  sites<-unique(days$site);hier<-split(seq_len(nrow(days)),days$site)
  hier<-lapply(hier,function(ix)split(ix,days$participant_key[ix]))
  W<-matrix(0,B,nrow(days));set.seed(cfg$seed+301L)
  for(b in seq_len(B)) {
    ix<-unlist(lapply(sample(sites,length(sites),replace=TRUE),function(site) {
      pp<-hier[[site]]
      unlist(lapply(sample.int(length(pp),length(pp),replace=TRUE),function(i)pp[[i]][sample.int(length(pp[[i]]),length(pp[[i]]),replace=TRUE)]),use.names=FALSE)
    }),use.names=FALSE)
    W[b,]<-tabulate(ix,nbins=nrow(days))
  }
  W
}
attribution_axis <- function(axis) {
  b<-readRDS(paste0('results/models/primary_',axis,'.rds'))
  # Sequential within an axis; each expensive subset has an atomic resume checkpoint.
  fits<-lapply(0:127,subset_fit,bundle=b)
  d<-b$data;days<-d |> distinct(site,participant_key,participant_day_key)
  dayid<-match(d$participant_day_key,days$participant_day_key)
  n<-tabulate(dayid,nbins=nrow(days));sy<-rowsum(d$response,dayid,reorder=FALSE)[,1]
  # rowsum output follows first appearance, the same order as distinct(days).
  sy2<-rowsum(d$response^2,dayid,reorder=FALSE)[,1]
  E<-vapply(fits,function(z)rowsum((d$response-z$prediction)^2,dayid,reorder=FALSE)[,1],numeric(nrow(days)))
  weights<-hierarchical_day_weights(days,cfg$bootstrap)
  N<-drop(weights%*%n);Y<-drop(weights%*%sy);Y2<-drop(weights%*%sy2);SST<-Y2-Y^2/N
  bootstrap_R2<-1-sweep(weights%*%E,1,SST,'/')
  M<-shapley_map();R2<-vapply(fits,`[[`,numeric(1),'R2');phi<-drop(R2%*%M);boot<-bootstrap_R2%*%M
  stopifnot(abs(sum(phi)-(R2[128]-R2[1]))<1e-9,
    max(abs(rowSums(boot)-(bootstrap_R2[,128]-bootstrap_R2[,1])))<1e-8)
  ci<-apply(boot,2,quantile,c(.025,.975),type=8)
  out<-tibble(axis=axis,component=blocks,allocation=phi,lower=ci[1,],upper=ci[2,],
    baseline_R2=R2[1],full_R2=R2[128],bootstrap_replicates=cfg$bootstrap,
    interval='Site-participant-day bootstrap of fixed predictions; conditional on fitted models')
  out<-bind_cols(out,count_sample(d)[rep(1,7),])
  save_table(out,paste0('shapley_',axis))
  design<-bind_rows(lapply(fits,function(z)tibble(axis=axis,mask=z$mask,included=paste(z$included,collapse=';'),formula=z$formula,
    R2=z$R2,rho=z$rho,converged=if(is.null(z$converged))NA else z$converged,warnings=paste(z$warnings,collapse='; '))))
  save_table(bind_cols(design,count_sample(d)[rep(1,128),]),paste0('shapley_subsets_',axis))
  atomic_rds(list(table=out,bootstrap=boot,baseline_bootstrap=bootstrap_R2[,1],full_bootstrap=bootstrap_R2[,128]),paste0('results/subsets/',axis,'_allocation.rds'))
  log_step('Completed 128 subsets and bootstrap:',axis)
  out
}
run_attribution <- function(prepared) {
  for(axis in c('clock','phase')) attribution_axis(axis)
}
