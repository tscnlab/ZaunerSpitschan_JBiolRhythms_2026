# Emens et al. (2009), doi:10.1080/07420520902821077.
# An exploratory balance of fitted exposure in two equal 12-hour regions.
emens_balance <- function(bundle, inference, draws=10000L, seed=cfg$seed+909L) {
  stopifnot(bundle$axis=='phase', draws>=1000L)
  x<-inference$cell_design$x
  X<-inference$pattern_matrices[grepl('^Q[14]_',names(inference$pattern_matrices))]
  stopifnot(length(X)==6L, length(x)==144L,
    all(abs(x-seq(-12,12-1/6,by=1/6))<1e-10))
  keep<-which(Reduce(`+`,lapply(X,function(m)colSums(abs(m))))>1e-12)
  beta<-coef(bundle$fit)[keep];cv<-covariance_info(bundle$fit)
  V<-cv$V[keep,keep,drop=FALSE]
  set.seed(seed)
  beta_draws<-mgcv::rmvn(draws,beta,V)
  eta<-lapply(X,function(m)drop(m[,keep,drop=FALSE]%*%beta))
  eta_draws<-lapply(X,function(m)m[,keep,drop=FALSE]%*%t(beta_draws))
  # This transforms the existing fitted curves, not the original observations
  # or a newly fitted log10(EDI+1) response. Put the matrix first in pmax so its
  # dimensions are retained. Negative back-transformed values are set to zero.
  emens_scale<-function(v)log10(pmax(10^v-.1,0)+1)
  y<-lapply(eta,emens_scale);ys<-lapply(eta_draws,emens_scale)
  counts<-inference$curves |> filter(quartile %in% c('Q1','Q4')) |>
    distinct(quartile,context,direction,participants,participant_days,half_hours) |>
    mutate(key=paste(quartile,context,sep='_'))
  summaries<-contrasts<-masks<-components<-checks<-list()
  add_check<-function(name,passed,value=NA_real_)checks[[length(checks)+1L]]<<-
    tibble(check=name,passed=isTRUE(passed),value=value)
  ci<-function(v)unname(quantile(v,c(.025,.975),type=8))
  for(offset in c(-6.95,-6.45,-5.95)) {
    phase<-(x-offset)%%24
    advance<-phase>=8 & phase<20 # DLMO -16 to -4 hours, modulo 24.
    delay<-!advance             # DLMO -4 to +8 hours.
    add_check(paste('Two equal 12-hour regions',offset),sum(advance)==72L && sum(delay)==72L)
    masks[[length(masks)+1L]]<-tibble(dlmo_offset=offset,x=x,
      phase_relative_to_dlmo=phase,region=ifelse(advance,'Advance','Delay'))
    point<-sim<-list()
    for(key in names(y)) {
      a<-mean(y[[key]][advance]);d<-mean(y[[key]][delay])
      ad<-colMeans(ys[[key]][advance,,drop=FALSE])
      dd<-colMeans(ys[[key]][delay,,drop=FALSE])
      stopifnot(d>0,all(is.finite(dd)),all(dd>0))
      point[[key]]<-a/d;sim[[key]]<-ad/dd;bounds<-ci(sim[[key]])
      row<-counts[counts$key==key,]
      summaries[[length(summaries)+1L]]<-bind_cols(row,tibble(dlmo_offset=offset,
        advance_mean_log10=a,delay_mean_log10=d,estimate=a/d,
        lower=bounds[1],upper=bounds[2],minimum_draw_denominator=min(dd),
        draws=draws,covariance=cv$label,
        interval_method='Joint coefficient draws; 95% percentile interval; alignment fixed'))
      # Equal-width four-hour bins reproduce the paper's three-bin means.
      bins<-floor(((phase+4)%%24)/4)
      means<-tapply(y[[key]],bins,mean)
      bin_ratio<-mean(means[as.character(3:5)])/mean(means[as.character(0:2)])
      add_check(paste('Six-bin and grid balances agree',offset,key),abs(bin_ratio-a/d)<1e-12,
        abs(bin_ratio-a/d))
      add_check(paste('Positive simulated delay denominator',offset,key),all(dd>0),min(dd))
      if(offset== -6.45)components[[length(components)+1L]]<-tibble(x=x,
        fitted_log10_edi_plus_01=eta[[key]],fitted_melanopic_edi=pmax(10^eta[[key]]-.1,0),
        emens_log10_edi_plus_1=y[[key]],region=ifelse(advance,'Advance','Delay'),key=key)
    }
    for(context in c('Overall','Work','Free')) {
      morning<-paste0('Q1_',context);evening<-paste0('Q4_',context)
      for(kind in c('Evening / morning','Evening minus morning')) {
        s<-if(kind=='Evening / morning')sim[[evening]]/sim[[morning]]else sim[[evening]]-sim[[morning]]
        value<-if(kind=='Evening / morning')point[[evening]]/point[[morning]]else point[[evening]]-point[[morning]]
        bounds<-ci(s)
        contrasts[[length(contrasts)+1L]]<-tibble(context=context,dlmo_offset=offset,
          contrast=kind,estimate=value,lower=bounds[1],upper=bounds[2],
          draws=draws,covariance=cv$label)
      }
    }
  }
  add_check('Zero light maps to zero',emens_scale(-1)==0)
  result<-list(summary=bind_rows(summaries),contrasts=bind_rows(contrasts),
    regions=bind_rows(masks),components=bind_rows(components),checks=bind_rows(checks),
    floor_counts=tibble(key=names(eta),negative_fitted_points=vapply(eta,function(v)sum(v< -1),integer(1)),
      total_points=length(x)),seed=seed,draws=draws)
  stopifnot(nrow(result$summary)==18L,nrow(result$contrasts)==18L,all(result$checks$passed))
  result
}

emens_balance_export <- function(result) {
  for(n in c('summary','contrasts','regions','components','checks','floor_counts'))
    save_table(result[[n]],paste0('emens_balance_',n),
      dir=if(n %in% c('checks','floor_counts'))'results/provenance'else'results/tables')
  atomic_rds(result,'results/data/emens_balance.rds')
  invisible(result)
}
