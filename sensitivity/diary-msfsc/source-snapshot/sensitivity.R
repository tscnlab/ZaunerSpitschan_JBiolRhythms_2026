sensitivity_predictions <- function(bundle,ref,free_weight) {
  c<-standardized_cells(bundle,ref);beta<-coef(bundle$fit);cr<-contrast_registry(free_weight)
  cells<-lapply(c$X,function(X)drop(X%*%beta));out<-list()
  for(i in seq_len(nrow(cr$W)))out[[i]]<-tibble(contrast=rownames(cr$W)[i],x=c$x,
    estimate=Reduce(`+`,Map(function(y,w)y*w,cells,cr$W[i,])))
  bind_rows(out)
}
run_sensitivity <- function(prepared) {
  summaries<-curves<-list()
  for(axis in c('clock','phase')) {
    b<-readRDS(paste0('results/models/primary_',axis,'.rds'));d<-b$data
    ref<-reference_days(d);free_weight<-sum(ref$w[ref$daytype=='Free'])
    original<-sensitivity_predictions(b,ref,free_weight)
    # Target the largest site-by-quartile imbalance and the largest participant residual RMSE.
    ss<-d |> distinct(site,participant_key,quartile) |> count(site,quartile) |> group_by(site) |>
      summarise(concentration=max(n)/sum(n),participants=sum(n)) |> arrange(desc(concentration),desc(participants))
    delete_site<-as.character(ss$site[1])
    pp<-d |> mutate(e=residuals(b$fit)) |> summarise(RMSE=sqrt(mean(e^2)),.by=participant_key) |> arrange(desc(RMSE))
    delete_person<-pp$participant_key[1]
    specifications<-list(
      larger_basis=list(raw=prepared$frames$glasses,k=18L,kp=14L,rho=b$rho),
      lower_rho=list(raw=prepared$frames$glasses,k=cfg$k,kp=cfg$kp,rho=max(-.95,b$rho-.1)),
      higher_rho=list(raw=prepared$frames$glasses,k=cfg$k,kp=cfg$kp,rho=min(.95,b$rho+.1)),
      delete_site=list(raw=filter(prepared$frames$glasses,as.character(site)!=delete_site),k=cfg$k,kp=cfg$kp,rho=b$rho),
      delete_participant=list(raw=filter(prepared$frames$glasses,participant_key!=delete_person),k=cfg$k,kp=cfg$kp,rho=b$rho),
      variance_weighted=list(raw=prepared$frames$glasses,k=cfg$k,kp=cfg$kp,rho=b$rho))
    for(name in names(specifications)) {
      spec<-specifications[[name]];path<-paste0('results/diagnostics/sensitivity_',axis,'_',name,'.rds')
      signature<-digest::digest(list(b$signature,spec,name,'sensitivity-v1'))
      old<-if(file.exists(path))readRDS(path)else NULL
      if(!is.null(old)&&identical(old$signature,signature))out<-old else {
        prep<-fit_preprocess(spec$raw,prepared$cutoffs);sd<-apply_preprocess(spec$raw,prep,axis)
        weights<-NULL;variance_model<-NULL
        if(name=='variance_weighted') {
          vd<-data.frame(mu=fitted(b$fit),e2=pmax(b$fit$std.rsd^2,1e-8))
          variance_model<-gam(e2~s(mu,bs='cr',k=8),data=vd,family=Gamma(link='log'),method='REML')
          weights<-1/as.numeric(predict(variance_model,type='response'));weights<-weights/mean(weights)
        }
        f<-fit_bam(sd,model_formula(k=spec$k,kp=spec$kp),spec$rho,paste('Sensitivity',axis,name),weights=weights)
        sb<-list(fit=f,data=sd,axis=axis,label=paste0('sensitivity_',name),rho=spec$rho,prep=prep)
        predictions<-sensitivity_predictions(sb,ref,free_weight)
        delta<-predictions |> left_join(rename(original,original=estimate),by=c('contrast','x'),relationship='one-to-one') |> mutate(change=estimate-original)
        stats<-bind_cols(tibble(axis=axis,sensitivity=name,maximum_contrast_change=max(abs(delta$change)),rms_contrast_change=sqrt(mean(delta$change^2)),
          rho=spec$rho,k=spec$k,participant_k=spec$kp,deleted=if(name=='delete_site')delete_site else if(name=='delete_participant')delete_person else'',
          convergence=isTRUE(f$converged),warnings=paste(attr(f,'analysis_warnings'),collapse='; ')),count_sample(sd))
        if(name=='variance_weighted') {
          # Check whether variance weighting reduces spread in the AR-standardized innovations.
          spread<-tibble(decile=dplyr::ntile(fitted(f),10),innovation=f$std.rsd) |> summarise(sd=sd(innovation),half_hours=n(),.by=decile)
          stats$innovation_spread_ratio<-max(spread$sd)/min(spread$sd)
          save_table(spread,paste0('variance_weighted_',axis,'_spread'),'results/diagnostics')
          atomic_rds(list(bundle=sb,variance_model=variance_model,weights=weights),paste0('results/models/variance_weighted_',axis,'.rds'))
        }
        out<-list(stats=stats,predictions=delta,signature=signature)
        atomic_rds(out,path);rm(f,sb);gc()
      }
      summaries[[paste(axis,name)]]<-out$stats
      curves[[paste(axis,name)]]<-mutate(out$predictions,axis=axis,sensitivity=name)
    }
    rm(b);gc()
  }
  save_table(bind_rows(summaries),'sensitivity_summary');save_table(bind_rows(curves),'sensitivity_contrasts')
  invisible(bind_rows(summaries))
}
