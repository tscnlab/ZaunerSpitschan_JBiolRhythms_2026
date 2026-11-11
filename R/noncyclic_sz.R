# Controlled sensitivity: retain the common cc curve and all fs/ti terms;
# replace only the four sz terms' cyclic cubic margins by ordinary cubic
# regression margins. This local constructor supplies their own full knot
# sequence because the shared cc boundary knots are not valid cr knots.
smooth.info.szcr.smooth.spec <- function(object) {
  object$tensor.possible <- TRUE
  object
}
smooth.construct.szcr.smooth.spec <- function(object, data, knots) {
  local_knots <- knots
  local_knots[['t']] <- seq(0, 24, length.out = object$bs.dim)
  mgcv:::smooth.construct.sz.smooth.spec(object, data, local_knots)
}

noncyclic_sz_formula <- function() {
  response ~ s(t, bs='cc', k=12) +
    s(t, quartile, bs='szcr', xt=list(bs='cr'), k=12) +
    s(t, daytype, bs='szcr', xt=list(bs='cr'), k=12) +
    s(t, quartile, daytype, bs='szcr', xt=list(bs='cr'), k=12) +
    s(t, site, bs='szcr', xt=list(bs='cr'), k=12) +
    s(photo, bs='cr', k=5) +
    ti(t, photo, bs=c('cc','cr'), k=c(12,5)) +
    s(t, participant, bs='fs', xt=list(bs='cc'), k=10) +
    s(participant_day, bs='re')
}

noncyclic_sz_fit <- function(axis) {
  baseline_path <- paste0('results/models/primary_', axis, '.rds')
  baseline <- readRDS(baseline_path)
  d <- baseline$data
  formula <- noncyclic_sz_formula()
  signature <- digest::digest(list(baseline=tools::md5sum(baseline_path),
    formula=deparse(formula), knots=seq(0,24,length.out=12), version=1L,
    mgcv=as.character(packageVersion('mgcv'))))
  path <- paste0('results/models/noncyclic_sz_',axis,'.rds')
  if(file.exists(path)) {
    cached <- readRDS(path)
    if(identical(cached$signature,signature)) return(cached)
  }
  preliminary <- fit_bam(d,formula,0,paste('noncyclic sz',axis,'preliminary'))
  set.seed(cfg$seed)
  preliminary_k <- mgcv::k.check(preliminary,subsample=5000,n.rep=100)
  preliminary_patterns <- tibble(t=d$t,residual=residuals(preliminary),
    quartile=d$quartile,daytype=d$daytype) |>
    mutate(hour=floor(t)) |>
    summarise(mean_residual=mean(residual),sd_residual=sd(residual),
      half_hours=n(),.by=c(hour,quartile,daytype))
  save_table(cbind(term=rownames(preliminary_k),as.data.frame(preliminary_k)),
    paste0('noncyclic_sz_',axis,'_preliminary_k'),'results/diagnostics')
  save_table(preliminary_patterns,paste0('noncyclic_sz_',axis,'_preliminary_patterns'),
    'results/diagnostics')
  rho <- pmin(.95,pmax(-.95,lag_cor(residuals(preliminary),d$AR_start)))
  rm(preliminary); gc()
  f <- fit_bam(d,formula,rho,paste('noncyclic sz',axis,'AR1'))
  out <- list(fit=f,data=d,prep=baseline$prep,axis=axis,label='noncyclic_sz',
    rho=rho,signature=signature,preliminary_k=preliminary_k)
  atomic_rds(out,path)
}

# Exact endpoint evaluations cannot use wrap24(24), which would hide a jump.
noncyclic_sz_cells_at <- function(bundle,t) {
  ref <- reference_days(bundle$data) |> summarise(w=sum(w),.by=photo)
  nd <- prediction_template(bundle,length(t)); nd$t <- t; nd$photo <- 0
  X0 <- population_matrix(bundle,nd); averaged <- X0
  for(sm in bundle$fit$smooth) if('photo' %in% sm$term) {
    ix <- sm$first.para:sm$last.para
    averaged[,ix] <- 0
    for(j in seq_len(nrow(ref))) {
      nd$photo <- ref$photo[j]
      averaged[,ix] <- averaged[,ix]+ref$w[j]*mgcv::PredictMat(sm,nd)
    }
  }
  nd$photo <- 0; delta <- averaged-X0; cells <- cell_registry()
  X <- lapply(seq_len(nrow(cells)),function(j) {
    z <- nd
    for(v in c('quartile','direction','extremity','daytype')) z[[v]] <- cells[[v]][j]
    population_matrix(bundle,z)+delta
  })
  names(X) <- cells$cell
  X
}

noncyclic_sz_boundary <- function(bundle) {
  X <- noncyclic_sz_cells_at(bundle,c(0,24))
  ref <- reference_days(bundle$data)
  fw <- sum(ref$w[ref$daytype=='Free'])
  cr <- contrast_registry(fw); V <- covariance_info(bundle$fit)$V
  beta <- coef(bundle$fit); curves <- contrasts <- list()
  for(context in c('Overall','Work','Free')) for(q in paste0('Q',1:4)) {
    cells <- if(context=='Overall')paste(q,c('Work','Free'),sep='_') else paste(q,context,sep='_')
    M <- if(context=='Overall')(1-fw)*X[[cells[1]]]+fw*X[[cells[2]]] else X[[cells]]
    # Forward traversal of the daily boundary: value immediately after minus before.
    jump <- linear_summary(M[1,,drop=FALSE]-M[2,,drop=FALSE],beta,V)
    curves[[length(curves)+1L]] <- bind_cols(tibble(model=bundle$label,axis=bundle$axis,
      context=context,quartile=q,before=drop(M[2,,drop=FALSE]%*%beta),
      after=drop(M[1,,drop=FALSE]%*%beta)),jump,counts_for_cells(bundle$data,cells))
  }
  for(i in seq_len(nrow(cr$W))) {
    M <- Reduce(`+`,Map(function(x,w)x*w,X,cr$W[i,]))
    jump <- linear_summary(M[1,,drop=FALSE]-M[2,,drop=FALSE],beta,V)
    w <- cr$W[i,]; a <- counts_for_cells(bundle$data,names(w)[w>0]); b <- counts_for_cells(bundle$data,names(w)[w<0])
    contrasts[[i]] <- bind_cols(tibble(model=bundle$label,axis=bundle$axis),cr$table[i,],
      tibble(before=drop(M[2,,drop=FALSE]%*%beta),after=drop(M[1,,drop=FALSE]%*%beta)),jump) |>
      mutate(positive_participants=a$participants,positive_days=a$participant_days,
        negative_participants=b$participants,negative_days=b$participant_days,
        half_hours=nrow(bundle$data))
  }
  list(curves=contrast_display(bind_rows(curves)),contrasts=contrast_display(bind_rows(contrasts)))
}

# Intervals are formed on the native [0,24) coordinate. The two ends are
# deliberately not joined: this sensitivity no longer imposes continuity there.
noncyclic_sz_windows <- function(inference,bundle) {
  out <- list(); d <- bundle$data
  for(id in inference$registry$table$contrast) {
    z <- filter(inference$contrasts,contrast==id) |> mutate(t=wrap24(x)) |> arrange(t)
    signs <- ifelse(z$significant,ifelse(z$estimate>0,1L,-1L),0L)
    rr <- rle(signs); ends <- cumsum(rr$lengths); starts <- ends-rr$lengths+1L
    w <- inference$registry$W[id,]
    for(j in which(rr$values!=0)) {
      ix <- seq.int(starts[j],ends[j])
      lo <- max(0,z$t[ix[1]]-cfg$step/2)
      hi <- min(24,z$t[tail(ix,1)]+cfg$step/2)
      inside <- d$t>=lo & d$t<hi
      a <- counts_for_cells(d[inside,],names(w)[w>0]); b <- counts_for_cells(d[inside,],names(w)[w<0])
      if(min(a$participants,b$participants)<cfg$support_n) next
      out[[length(out)+1L]] <- tibble(model=bundle$label,axis=bundle$axis,contrast=id,
        description=z$description[1],context=z$context[1],direction=if(rr$values[j]>0)'Positive'else'Negative',
        native_start=lo,native_end=hi,start=if(bundle$axis=='phase')wrap12(lo)else lo,
        end=if(bundle$axis=='phase')wrap12(hi)else hi,duration_hours=hi-lo,
        mean_difference=mean(z$estimate[ix]),mean_ratio=10^mean(z$estimate[ix]),
        minimum_difference=min(z$estimate[ix]),maximum_difference=max(z$estimate[ix]),
        positive_participants=a$participants,positive_days=a$participant_days,positive_half_hours=a$half_hours,
        negative_participants=b$participants,negative_days=b$participant_days,negative_half_hours=b$half_hours)
    }
  }
  if(length(out))bind_rows(out)else tibble(model=character(),axis=character(),contrast=character())
}

noncyclic_sz_acceptance <- function(bundle,baseline) {
  checks <- list()
  add <- function(name,error,tol=1e-8) checks[[length(checks)+1L]] <<-
    tibble(axis=bundle$axis,check=name,error=as.numeric(error),tolerance=tol,passed=is.finite(error)&&error<=tol)
  add('Identical analysis data and row ordering',as.numeric(!identical(bundle$data,baseline$data)),0)
  add('Identical preprocessing and cutoffs',as.numeric(!identical(bundle$prep,baseline$prep)),0)
  add('Sample: 140 participants, 778 days, 35949 half-hours',
    sum(abs(unlist(count_sample(bundle$data))-c(140,778,35949))),0)
  labs <- vapply(bundle$fit$smooth,`[[`,character(1),'label')
  base_labs <- vapply(baseline$fit$smooth,`[[`,character(1),'label')
  add('Same smooth term labels',as.numeric(!identical(labs,base_labs)),0)
  nd <- prediction_template(bundle,7); nd$t <- c(0,1,6,12,18,23,24)
  nd$photo <- seq(min(bundle$data$photo),max(bundle$data$photo),length.out=7)
  for(i in seq_along(bundle$fit$smooth)) {
    sm <- bundle$fit$smooth[[i]]; old <- baseline$fit$smooth[[i]]
    if(!inherits(sm,'sz.interaction')) {
      # bam can rotate/center coefficient bases when another term changes.
      # Compare their defining basis classes, dimensions and knots, not the
      # arbitrary coefficient coordinates of those representations.
      definition <- function(s) list(class=class(s),term=s$term,bs.dim=s$bs.dim,
        xp=s$xp,base=s$base$bs,margin=lapply(s$margin,function(m)
          list(class=class(m),term=m$term,bs.dim=m$bs.dim,xp=m$xp)))
      add(paste('Unchanged basis definition',sm$label),
        as.numeric(!isTRUE(all.equal(definition(sm),definition(old),tolerance=1e-12))),0)
      if('t' %in% sm$term) {
        endpoint <- prediction_template(bundle,2); endpoint$t <- c(0,24)
        endpoint$photo <- 0.5
        P <- mgcv::PredictMat(sm,endpoint)
        add(paste('Retained cyclic endpoint',sm$label),max(abs(P[1,]-P[2,])))
      }
    } else {
      add(paste('Noncyclic cr margin',sm$label),as.numeric(!identical(sm$base$bs,'cr.smooth')),0)
      for(v in intersect(c('quartile','daytype','site'),sm$term)) {
        summed <- Reduce(`+`,lapply(levels(bundle$data[[v]]),function(level) {
          zz <- nd; zz[[v]] <- factor(level,levels=levels(bundle$data[[v]]))
          P <- predict(bundle$fit,zz,type='lpmatrix',discrete=FALSE)
          P[,sm$first.para:sm$last.para,drop=FALSE]
        }))
        add(paste('Sum-to-zero',sm$label,'across',v),max(abs(summed)))
      }
    }
  }
  out <- bind_rows(checks)
  save_table(out,paste0('noncyclic_sz_',bundle$axis,'_acceptance'),'results/provenance')
  out
}

noncyclic_sz_pattern_plot <- function(curves,context) {
  panels <- list()
  for(axis in c('clock','phase')) {
    z <- filter(curves,.data$axis==.env$axis,.data$context==.env$context) |>
      mutate(extremity=factor(extremity,c('Moderate','Extreme')),
        segment=if_else(axis=='phase'&x<0,'Before MSFsc','After MSFsc'))
    labs <- distinct(z,extremity,direction,participants,participant_days) |>
      mutate(label=paste0(direction,': ',participants,' participants; ',participant_days,' participant-days'),
        v=if_else(direction=='Morning',1.2,2.7))
    p <- ggplot(z,aes(x,estimate_melEDI,color=direction,linetype=model,
      group=interaction(direction,model,segment)))
    if(axis=='phase')p <- p+dlmo_reference_band()+geom_vline(xintercept=0,color='grey45',linetype=3)
    p <- p+geom_ribbon(data=filter(z,model=='Noncyclic sz'),aes(ymin=lower_melEDI,ymax=upper_melEDI,fill=direction),
      alpha=.13,color=NA,linetype=0)+geom_line(linewidth=.65)+facet_wrap(~extremity,nrow=1)+
      geom_text(data=labs,aes(x=Inf,y=Inf,label=label,color=direction,vjust=v),inherit.aes=FALSE,hjust=1.03,size=3,show.legend=FALSE)+
      scale_color_manual(values=c(Morning='#0072B2',Evening='#D55E00'))+
      scale_fill_manual(values=c(Morning='#0072B2',Evening='#D55E00'),guide='none')+
      scale_linetype_manual(values=c('Cyclic sz'='dashed','Noncyclic sz'='solid'))+
      time_scale(axis)+intensity_scale(expand=expansion(mult=c(.03,.23)))+
      labs(x=axis_title(axis),y='Melanopic EDI (lx; symlog)',color='Chronotype',linetype='Specification',tag=if(axis=='clock')'A'else'B')
    panels[[axis]] <- p
  }
  patchwork::wrap_plots(panels,ncol=1,guides='collect')+
    patchwork::plot_annotation(title=paste('Near-eye:',context,'population curves'),
      subtitle='Matched observations and standardization; only sz time margins change',
      caption='A: clock time. B: MSFsc-relative time. Shading: 95% pointwise bands for noncyclic sz curves.\nCyclic curves are dashed. Phase lines are split at MSFsc (0 h), the unconstrained join. Striped region: published DLMO reference.')
}

noncyclic_sz_contrast_plot <- function(contrasts,context='Overall') {
  panels <- list()
  for(axis in c('clock','phase')) {
    z <- filter(contrasts,.data$axis==.env$axis,.data$context==.env$context,family=='Evening minus morning') |>
      mutate(description=factor(gsub(' minus ',' / ',description,fixed=TRUE),
        levels=c('Evening / morning: Moderate','Evening / morning: Extreme')),
        segment=if_else(axis=='phase'&x<0,'Before MSFsc','After MSFsc'))
    labs <- distinct(z,description,positive_participants,positive_days,negative_participants,negative_days) |>
      mutate(label=paste0('Evening: ',positive_participants,' participants; ',positive_days,' participant-days\nMorning: ',negative_participants,' participants; ',negative_days,' participant-days'))
    p <- ggplot(z,aes(x,ratio,linetype=model,group=interaction(model,segment)))
    if(axis=='phase')p <- p+dlmo_reference_band()+geom_vline(xintercept=0,color='grey45',linetype=3)
    p <- p+geom_hline(yintercept=1,color='grey55')+
      geom_ribbon(data=filter(z,model=='Noncyclic sz'),aes(ymin=ratio_lower,ymax=ratio_upper),
        fill='#0072B2',alpha=.13,linetype=0)+geom_line(color='#0072B2')+
      geom_point(data=filter(z,model=='Noncyclic sz',significant),color='#B32025',size=.65)+
      geom_text(data=labs,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=1.03,vjust=1.05,size=2.8)+
      facet_wrap(~description,nrow=1)+scale_linetype_manual(values=c('Cyclic sz'='dashed','Noncyclic sz'='solid'))+
      time_scale(axis)+factor_scale(expand=expansion(mult=c(.04,.25)))+
      labs(x=axis_title(axis),y='Evening / morning factor for (melEDI + 0.1)',linetype='Specification',tag=if(axis=='clock')'A'else'B')
    panels[[axis]] <- p
  }
  patchwork::wrap_plots(panels,ncol=1,guides='collect')+
    patchwork::plot_annotation(title=paste('Near-eye:',context,'evening / morning contrasts'),
      subtitle='Same 20 registered contrasts; 10,000 coefficient draws for the new specification',
      caption=paste0('A: clock time. B: MSFsc-relative time. Bands for noncyclic sz: simultaneous across 20 contrasts and the grid,\n95% within each time coordinate. Red points: candidate exclusions of 1 with at least ',cfg$support_n,' people per cell.\nThe Gaussian response has a floor at zero exposure. The phase join is at 0 h; lines and ribbons do not bridge it.'))
}

noncyclic_sz_outputs <- function() {
  models <- list(); inferences <- list(); verdicts <- list(); checks <- list()
  for(axis in c('clock','phase')) {
    b <- noncyclic_sz_fit(axis); base <- readRDS(paste0('results/models/primary_',axis,'.rds'))
    checks[[axis]] <- noncyclic_sz_acceptance(b,base)
    diagnostic_cache <- paste0('results/diagnostics/noncyclic_sz_',axis,'_diagnostics.rds')
    if(file.exists(diagnostic_cache)) dg <- readRDS(diagnostic_cache) else {
      dg <- diagnose_model(b); atomic_rds(dg,diagnostic_cache)
    }
    inference_cache <- paste0('results/models/noncyclic_sz_',axis,'_inference.rds')
    inf <- if(file.exists(inference_cache)) readRDS(inference_cache) else NULL
    if(is.null(inf) || !identical(inf$axis_coverage,cfg$curve_coverage) ||
       !identical(inf$minimum_support,cfg$support_n)) {
      inf <- infer_curves(b)
      inf$windows <- noncyclic_sz_windows(inf,b)
      save_table(inf$windows,paste0('noncyclic_sz_',axis,'_significant_windows'))
      atomic_rds(inf,inference_cache)
    }
    models[[paste0(axis,'_base')]] <- base; models[[axis]] <- b
    inferences[[axis]] <- inf; verdicts[[axis]] <- dg$flags
  }
  checks <- bind_rows(checks); stopifnot(all(checks$passed))
  fit_table <- bind_rows(lapply(models,model_summary))
  save_table(fit_table,'noncyclic_sz_fit_comparison')
  baseline_verdicts <- bind_rows(lapply(c('clock','phase'),function(axis)
    read.csv(paste0('results/diagnostics/primary_',axis,'_verdict.csv'))))
  diagnostic_comparison <- bind_rows(baseline_verdicts,bind_rows(verdicts))
  stopifnot(nrow(diagnostic_comparison)==4L,
    identical(names(diagnostic_comparison),names(baseline_verdicts)),
    !anyDuplicated(diagnostic_comparison[c('model','axis')]),
    !anyNA(diagnostic_comparison),
    all(vapply(diagnostic_comparison[vapply(diagnostic_comparison,is.numeric,logical(1))],
      function(x)all(is.finite(x)),logical(1))))
  save_table(diagnostic_comparison,'noncyclic_sz_diagnostic_comparison')
  boundaries <- lapply(models,noncyclic_sz_boundary)
  boundary_curves <- bind_rows(lapply(boundaries,`[[`,'curves'))
  boundary_contrasts <- bind_rows(lapply(boundaries,`[[`,'contrasts'))
  save_table(boundary_curves,'noncyclic_sz_boundary_curves')
  save_table(boundary_contrasts,'noncyclic_sz_boundary_contrasts')
  allcurves <- allcontrasts <- list()
  for(axis in c('clock','phase')) {
    allcurves[[axis]] <- bind_rows(
      read.csv(paste0('results/tables/primary_',axis,'_curves.csv')) |> mutate(model='Cyclic sz'),
      intensity_display(inferences[[axis]]$curves) |> mutate(model='Noncyclic sz')) |> mutate(axis=axis)
    allcontrasts[[axis]] <- bind_rows(
      read.csv(paste0('results/tables/primary_',axis,'_contrasts.csv')) |> mutate(model='Cyclic sz'),
      inferences[[axis]]$contrasts |> mutate(model='Noncyclic sz')) |> mutate(axis=axis)
  }
  allcurves <- bind_rows(allcurves); allcontrasts <- bind_rows(allcontrasts)
  save_table(allcurves,'noncyclic_sz_curve_comparison')
  save_table(allcontrasts,'noncyclic_sz_contrast_comparison')
  for(context in c('Overall','Work','Free')) {
    save_plot(noncyclic_sz_pattern_plot(allcurves,context),paste0('noncyclic_sz_patterns_',tolower(context)),12,9)
    save_plot(noncyclic_sz_contrast_plot(allcontrasts,context),paste0('noncyclic_sz_contrasts_',tolower(context)),12,9)
  }
  # CSV decimal rounding can differ from the newly generated floating-point
  # grid. Join by the exact integer 10-minute grid index, not decimal hours.
  paired <- allcontrasts |>
    mutate(grid_index=as.integer(round(x/cfg$step))) |>
    select(axis,contrast,grid_index,model,estimate,significant,supported)
  stopifnot(!anyDuplicated(paired[c('axis','contrast','grid_index','model')]))
  paired <- pivot_wider(paired,names_from=model,values_from=c(estimate,significant,supported))
  stopifnot(nrow(paired)==2L*20L*144L,!anyNA(paired),
    all(is.finite(paired$`estimate_Cyclic sz`)),all(is.finite(paired$`estimate_Noncyclic sz`)))
  stability <- paired |>
    summarise(maximum_abs_difference=max(abs(`estimate_Noncyclic sz`-`estimate_Cyclic sz`)),
      median_abs_difference=median(abs(`estimate_Noncyclic sz`-`estimate_Cyclic sz`)),
      cyclic_candidate_points=sum(`significant_Cyclic sz`),noncyclic_candidate_points=sum(`significant_Noncyclic sz`),
      candidate_agreement=mean(`significant_Cyclic sz`==`significant_Noncyclic sz`),
      participants=140L,participant_days=778L,half_hours=35949L,.by=c(axis,contrast))
  stopifnot(nrow(stability)==40L,all(vapply(stability[3:7],function(x)all(is.finite(x)),logical(1))))
  save_table(tibble(check=c('One observation per model, contrast and grid index',
    '5760 complete paired grid rows','40 finite contrast summary rows'),passed=TRUE,
    participants=140L,participant_days=778L,half_hours=35949L),
    'noncyclic_sz_comparison_acceptance','results/provenance')
  save_table(stability,'noncyclic_sz_contrast_stability')
  inputs <- c('results/data/prepared.rds','results/models/primary_clock.rds','results/models/primary_phase.rds',
    'results/models/noncyclic_sz_clock.rds','results/models/noncyclic_sz_phase.rds','R/noncyclic_sz.R','R/run_noncyclic_sz.R',
    'noncyclic-sz-sensitivity.qmd')
  save_table(tibble(path=inputs,md5=unname(tools::md5sum(inputs))),'noncyclic_sz_inputs','results/provenance')
  capture.output({
    cat('Commands:\nRscript --vanilla R/run_noncyclic_sz.R\nquarto render noncyclic-sz-sensitivity.qmd\n\n')
    print(Sys.time()); print(sessionInfo()); print(noncyclic_sz_formula())
    cat('\nCustom sz knots:',seq(0,24,length.out=12),'\n')
    cat('OpenMP:',Sys.getenv('MELIDOS_OPENMP_LIBRARY','/opt/homebrew/opt/libomp/lib/libomp.dylib'),'\n')
    print(fit_table); print(checks)
  },file='results/provenance/noncyclic_sz_runtime.txt')
  invisible(list(fit=fit_table,checks=checks,boundaries=boundary_curves))
}
