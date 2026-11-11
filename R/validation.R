validation_splits <- function(d) {
  set.seed(cfg$seed+101L)
  participants<-d |> distinct(site,participant_key) |> arrange(site,participant_key) |>
    group_by(site) |> mutate(fold=sample(rep(1:5,length.out=n()))) |> ungroup()
  splits<-lapply(1:5,function(i)list(id=paste0('new_participant_fold_',i),target='New participants',
    assessment=d$participant_key %in% participants$participant_key[participants$fold==i]))
  final<-d |> distinct(participant_key,participant_day_key,local_date) |> group_by(participant_key) |>
    filter(n()>=3,local_date==max(local_date)) |> ungroup()
  splits[[6]]<-list(id='known_participant_final_day',target='Known participants: final day',assessment=d$participant_day_key %in% final$participant_day_key)
  save_table(participants,'participant_folds','results/validation')
  save_table(final,'held_out_final_days','results/validation')
  splits
}
predict_assessment <- function(bundle,test,target) {
  nd<-apply_preprocess(test,bundle$prep,bundle$axis)
  exclude<-exclude_labels(bundle$fit,site=FALSE,participant=target=='New participants',day=TRUE)
  new_p<-!as.character(nd$participant) %in% levels(bundle$data$participant)
  if(any(new_p))nd$participant<-as.character(nd$participant)
  if(any(new_p))nd$participant[new_p]<-levels(bundle$data$participant)[1]
  nd$participant_day<-levels(bundle$data$participant_day)[1]
  stopifnot(all(as.character(nd$site) %in% levels(bundle$data$site)))
  nd<-set_factors(nd,bundle)
  pred<-as.numeric(predict(bundle$fit,newdata=nd,type='response',exclude=exclude,discrete=FALSE))
  stopifnot(all(is.finite(pred)))
  tibble(participant_key=nd$participant_key,participant_day_key=nd$participant_day_key,
    site=as.character(nd$site),local_date=nd$local_date,clock_bin=nd$clock_bin,observed=nd$response,predicted=pred,
    baseline=mean(bundle$data |> summarise(m=mean(response),.by=participant_key) |> pull(m)),
    photo_outside_training=nd$photo<min(bundle$data$photo)|nd$photo>max(bundle$data$photo))
}
score_rows <- function(d) {
  pp<-d |> summarise(mse=mean((predicted-observed)^2),mae=mean(abs(predicted-observed)),bias=mean(predicted-observed),
    null_mse=mean((baseline-observed)^2),.by=participant_key)
  w<-1/as.numeric(table(d$participant_key)[d$participant_key]); w<-w/sum(w)
  cal<-lm(observed~predicted,data=d,weights=w)
  bind_cols(tibble(RMSE=sqrt(mean(pp$mse)),MAE=mean(pp$mae),bias=mean(pp$bias),predictive_R2=1-mean(pp$mse)/mean(pp$null_mse),
    calibration_intercept=coef(cal)[1],calibration_slope=coef(cal)[2]),count_sample(d))
}
paired_validation <- function(predictions) {
  z<-predictions |> summarise(mse=mean((predicted-observed)^2),mae=mean(abs(predicted-observed)),bias=mean(predicted-observed),
    null_mse=mean((baseline-observed)^2),.by=c(target,axis,site,participant_key))
  out<-list();set.seed(cfg$seed+102L)
  for(target in unique(z$target)) {
    a<-filter(z,.data$target==.env$target,axis=='clock') |> arrange(site,participant_key)
    b<-filter(z,.data$target==.env$target,axis=='phase') |> arrange(site,participant_key)
    stopifnot(identical(a$participant_key,b$participant_key),isTRUE(all.equal(a$null_mse,b$null_mse)))
    score<-function(x,ix)c(RMSE=sqrt(mean(x$mse[ix])),MAE=mean(x$mae[ix]),bias=mean(x$bias[ix]),predictive_R2=1-mean(x$mse[ix])/mean(x$null_mse[ix]))
    delta<-score(b,seq_len(nrow(b)))-score(a,seq_len(nrow(a)))
    strata<-split(seq_len(nrow(a)),a$site)
    draws<-replicate(cfg$bootstrap,{ix<-unlist(lapply(strata,function(i)i[sample.int(length(i),length(i),replace=TRUE)]));score(b,ix)-score(a,ix)})
    ci<-apply(draws,1,quantile,c(.025,.975),type=8)
    out[[target]]<-bind_cols(tibble(target=target,metric=names(delta),phase_minus_clock=as.numeric(delta),lower=ci[1,],upper=ci[2,],
      bootstrap_replicates=cfg$bootstrap),count_sample(filter(predictions,.data$target==.env$target,axis=='clock'))[rep(1,4),])
  }
  bind_rows(out)
}
run_validation <- function(prepared) {
  d<-prepared$frames$glasses;splits<-validation_splits(d);preds<-counts<-list()
  for(split in splits) {
    train<-d[!split$assessment,];test<-d[split$assessment,]
    stopifnot(length(intersect(train$participant_day_key,test$participant_day_key))==0)
    if(split$target=='New participants')stopifnot(length(intersect(train$participant_key,test$participant_key))==0)
    prep<-fit_preprocess(train)
    for(axis in c('clock','phase')) {
      path<-paste0('results/validation/',split$id,'_',axis,'.rds')
      signature<-digest::digest(list(train,test,prep,cfg$version,deparse(model_formula())))
      old<-if(file.exists(path))readRDS(path)else NULL
      if(!is.null(old)&&identical(old$signature,signature))out<-old else {
        b<-fit_pair(train,prep,axis,split$id)
        pr<-predict_assessment(b,test,split$target) |> mutate(target=split$target,split=split$id,axis=axis)
        out<-list(predictions=pr,signature=signature,training=count_sample(train),assessment=count_sample(test),
          cutoffs=prep$cutoffs,centers=prep$centers,rho=b$rho,diagnostics=model_summary(b))
        atomic_rds(out,path);rm(b);gc()
      }
      preds[[paste(split$id,axis)]]<-out$predictions
      counts[[paste(split$id,axis)]]<-tibble(target=split$target,axis=axis,split=split$id,
        training_participants=out$training$participants,training_days=out$training$participant_days,training_half_hours=out$training$half_hours,
        assessment_participants=out$assessment$participants,assessment_days=out$assessment$participant_days,assessment_half_hours=out$assessment$half_hours,
        cutoff25=out$cutoffs[1],cutoff50=out$cutoffs[2],cutoff75=out$cutoffs[3],rho=out$rho,
        photo_outside_training=sum(out$predictions$photo_outside_training))
    }
  }
  predictions<-bind_rows(preds);save_table(predictions,'held_out_predictions','results/validation')
  save_table(bind_rows(counts),'validation_samples')
  scores<-predictions |> group_by(target,axis) |> group_modify(~score_rows(.x)) |> ungroup()
  folds<-predictions |> group_by(target,split,axis) |> group_modify(~score_rows(.x)) |> ungroup()
  save_table(scores,'validation_scores');save_table(folds,'validation_fold_scores')
  save_table(paired_validation(predictions),'validation_paired_differences')
  log_step('Completed predictive validation')
  invisible(scores)
}
