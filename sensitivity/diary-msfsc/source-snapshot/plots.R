dlmo_reference <- function() {
  # Facer-Childs et al. (2019), Table 1: early and late group means.
  z<-tibble(group=c('Early','Late'),MSFsc=c(2+24/60,6+52/60),DLMO=c(20+27/60,23+55/60)) |>
    mutate(relative_DLMO=wrap12(DLMO-MSFsc))
  save_table(z,'published_dlmo_reference')
  range(z$relative_DLMO)
}
axis_title <- function(axis) if(axis=='clock')'Local clock time'else'Hours relative to MSFsc'
# Reader-facing labels leave model factors and contrast definitions unchanged.
chronotype_plot_label <- function(x) {
  x<-as.character(x)
  for(pair in list(
                  c('Direction','Chronotype direction'),c('direction','chronotype direction'),
                  c('Extremity','Chronotype strength'),c('extremity','chronotype strength')))
    x<-gsub(paste0('\\b',pair[1],'\\b'),pair[2],x,perl=TRUE)
  gsub('Chronotype chronotype','Chronotype',x,fixed=TRUE)
}

# Grid hatching uses physical units, so diagonal stripes retain their angle and
# spacing in narrow facets and in both raster and vector output. The nested
# viewport clips stripes to the reference interval. Annotation grobs do not
# train either scale, which preserves log-factor and symlog-intensity limits.
makeContent.dlmo_reference_grob <- function(x) {
  width<-grid::convertWidth(grid::unit(1,'npc'),'mm',valueOnly=TRUE)
  height<-grid::convertHeight(grid::unit(1,'npc'),'mm',valueOnly=TRUE)
  intercept<-seq(-height,width,by=2)
  children<-grid::gTree(vp=grid::viewport(clip='on'),children=grid::gList(
    grid::rectGrob(gp=grid::gpar(fill=grDevices::adjustcolor('#D7C899',alpha.f=.68),col=NA)),
    grid::segmentsGrob(x0=grid::unit(intercept,'mm'),y0=grid::unit(0,'mm'),
      x1=grid::unit(intercept+height,'mm'),y1=grid::unit(height,'mm'),
      gp=grid::gpar(col=grDevices::adjustcolor('#8B7541',alpha.f=.72),lwd=.55))))
  grid::setChildren(x,grid::gList(children))
}
dlmo_reference_band <- function() {
  # Facer-Childs et al. (2019), Table 1; identical endpoints to dlmo_reference().
  limits<-range(wrap12(c(20+27/60,23+55/60)-c(2+24/60,6+52/60)))
  band_geom<-ggproto('GeomDlmoReference',Geom,
    required_aes=character(0),draw_key=draw_key_blank,
    draw_panel=function(data,panel_params,coord,xmin,xmax) {
      stopifnot(inherits(coord,'CoordCartesian'))
      # Use the panel viewport for vertical extent, rather than sending -Inf
      # through a logarithmic scale. The x coordinate is linear in phase hours.
      edges<-coord$transform(data.frame(x=c(xmin,xmax),y=c(0,0)),panel_params)$x
      grid::gTree(cl='dlmo_reference_grob',
        vp=grid::viewport(x=mean(edges),y=.5,width=diff(edges),height=1))
    })
  layer(geom=band_geom,stat='identity',position='identity',
    data=data.frame(.band=1),mapping=NULL,inherit.aes=FALSE,show.legend=FALSE,
    params=list(xmin=limits[1],xmax=limits[2]))
}
primary_phase_interval_band <- function(boundaries=TRUE) {
  selected<-read.csv('results/tables/primary_phase_significant_windows.csv') |>
    filter(contrast=='evening_morning_Extreme_Overall') |> arrange(start)
  stopifnot(nrow(selected)==2L)
  # Panel coordinates preserve the vertical extent on linear, log and symlog
  # axes without changing their limits. These are the primary-model intervals,
  # not intervals selected from the exploratory photoperiod comparison.
  band_geom<-ggproto('GeomPrimaryPhaseIntervals',Geom,
    required_aes=character(0),draw_key=draw_key_blank,
    draw_panel=function(data,panel_params,coord,starts,ends,boundaries) {
      stopifnot(inherits(coord,'CoordCartesian'))
      left<-coord$transform(data.frame(x=starts,y=0),panel_params)$x
      right<-coord$transform(data.frame(x=ends,y=0),panel_params)$x
      grid::grobTree(
        grid::rectGrob(x=(left+right)/2,y=.5,width=right-left,height=1,
          gp=grid::gpar(fill=grDevices::adjustcolor('grey65',alpha.f=.2),col=NA)),
        if(boundaries) grid::segmentsGrob(x0=c(left,right),x1=c(left,right),y0=0,y1=1,
          gp=grid::gpar(col='grey45',lty='dashed',lwd=1)) else grid::nullGrob())
    })
  layer(geom=band_geom,stat='identity',position='identity',
    data=data.frame(.band=1),mapping=NULL,inherit.aes=FALSE,show.legend=FALSE,
    params=list(starts=selected$start,ends=selected$end,boundaries=boundaries))
}

time_scale <- function(axis) if(axis=='clock') scale_x_continuous(breaks=seq(0,24,4),labels=function(x)ifelse(x==24,'24:00',clock_label(x)),limits=c(0,24),oob=scales::oob_keep,expand=expansion(mult=.025))else scale_x_continuous(breaks=seq(-12,12,4),limits=c(-12,12),oob=scales::oob_keep,expand=expansion(mult=.025))
intensity_scale <- function(...) scale_y_continuous(
  transform=LightLogR::symlog_trans(base=10,thr=1,scale=1),
  breaks=c(0,1,10,100,1000,10000,100000),labels=scales::label_number(big.mark=','),...)
factor_scale <- function(...) scale_y_log10(
  breaks=function(limits) if(diff(log10(limits))<=2)
    c(.01,.02,.03,.05,.1,.2,.3,.5,.75,1,1.5,2,3,5,10,20,30,50,100)else
    c(.001,.01,.1,.5,1,2,10,100,1000,10000),
  labels=function(x)vapply(x,function(v)format(v,trim=TRUE,scientific=FALSE,big.mark=','),character(1)),...)
pattern_figure <- function(curves,context,placement) {
  curves<-intensity_display(curves)
  band_description<-if('band_type'%in%names(curves)&&all(curves$band_type=='simultaneous'))
    '95% simultaneous bands across quartiles, day contexts and times within each coordinate'else'95% pointwise bands'
  panels<-list()
  for(axis in c('clock','phase')) {
    z<-filter(curves,.data$context==.env$context,.data$axis==.env$axis) |> mutate(extremity=factor(extremity,c('Moderate','Extreme')))
    stopifnot(n_distinct(z$axis)==1L,n_distinct(z$context)==1L,nrow(z)==4*length(unique(z$x)))
    lab<-distinct(z,extremity,direction,participants,participant_days) |> mutate(label=paste0(chronotype_plot_label(direction),': ',participants,' participants; ',participant_days,' participant-days'),v=if_else(direction=='Morning',1.2,2.7))
    p<-ggplot(z,aes(x,estimate_melEDI,color=direction,fill=direction))
    if(axis=='phase')p<-p+dlmo_reference_band()
    p<-p+geom_ribbon(aes(ymin=lower_melEDI,ymax=upper_melEDI),alpha=.16,color=NA)+geom_line(linewidth=.7)+facet_wrap(~extremity,nrow=1,labeller=as_labeller(function(x)paste('Chronotype strength:',tolower(x))))+
      geom_text(data=lab,aes(x=Inf,y=Inf,label=label,color=direction,vjust=v),inherit.aes=FALSE,hjust=1.03,size=3.1)+
      scale_color_manual(values=c(Morning='#0072B2',Evening='#D55E00'))+scale_fill_manual(values=c(Morning='#0072B2',Evening='#D55E00'))+
      time_scale(axis)+intensity_scale(limits=range(c(0,curves$lower_melEDI,curves$upper_melEDI)),expand=expansion(mult=c(.03,.24)))+
      labs(x=axis_title(axis),y='Melanopic EDI (lx; symlog)',color=NULL,fill=NULL,tag=if(axis=='clock')'A'else'B')+
      theme(legend.position='none',plot.tag=element_text(face='bold'))
    panels<-c(panels,list(p))
  }
  patchwork::wrap_plots(panels,ncol=1)+
    patchwork::plot_annotation(title=paste(placement,context,'daily exposure patterns'),
      subtitle='Common photoperiod distribution; site, participant and day effects excluded',
      caption=paste0('A: local clock time. B: circadian phase. Chronotype direction: morning (blue) and evening (orange).\nBack-transformed log-scale means and ',band_description,'.\nSymlog: linear below 1 lx, base-10 above 1 lx. Striped band: DLMO reference (Facer-Childs et al., 2019).\nOverall uses the same work/free-day mixture for every group.'))
}
difference_figure <- function(contrasts,family,context=NULL,placement='Near-eye') {
  z<-contrast_display(filter(contrasts,.data$family==.env$family))
  if(!is.null(context))z<-filter(z,.data$context==.env$context)
  z$description<-chronotype_plot_label(gsub(' minus ',' / ',z$description,fixed=TRUE))
  z$description<-factor(z$description,levels=unique(z$description))
  panels<-list()
  factor_limits<-range(c(1,z$ratio_lower,z$ratio_upper))
  numerator<-switch(family,'Evening minus morning'='Evening','Extreme minus moderate'='Extreme','Free minus work'='Free','Interactions'='Numerator cells')
  denominator<-switch(family,'Evening minus morning'='Morning','Extreme minus moderate'='Moderate','Free minus work'='Work','Interactions'='Denominator cells')
  for(axis in c('clock','phase')) {
    a<-filter(z,.data$axis==.env$axis)
    stopifnot(n_distinct(a$axis)==1L,nrow(a)==n_distinct(a$description)*length(unique(a$x)))
    lab<-a |> distinct(description,positive_participants,positive_days,negative_participants,negative_days) |>
      mutate(label=paste0(numerator,': ',positive_participants,' participants, ',positive_days,' participant-days\n',
        denominator,': ',negative_participants,' participants, ',negative_days,' participant-days'))
    p<-ggplot(a,aes(x,ratio))
    if(axis=='phase')p<-p+dlmo_reference_band()
    p<-p+geom_hline(yintercept=1,color='grey50')+geom_ribbon(aes(ymin=ratio_lower,ymax=ratio_upper),fill='#0072B2',alpha=.15)+
      geom_line(color='#0072B2')+geom_point(data=filter(a,significant),color='#B32025',size=.7)+
      geom_text(data=lab,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=1.02,vjust=1.1,size=2.8)+
      facet_wrap(~description,nrow=if(n_distinct(a$description)>2)2 else 1)+time_scale(axis)+
      factor_scale(limits=factor_limits,oob=scales::oob_keep,expand=expansion(mult=c(.05,.18)))+
      labs(x=axis_title(axis),y=if(family=='Interactions')'Relative factor (ratio of ratios)'else'Factor for (melEDI + 0.1)',tag=if(axis=='clock')'A'else'B')+
      theme(plot.tag=element_text(face='bold'))
    panels<-c(panels,list(p))
  }
  reference<-switch(family,'Evening minus morning'='Reference: morning group','Extreme minus moderate'='Reference: moderate group',
    'Free minus work'='Reference: work days','Interactions'='Reference: no interaction (relative ratio = 1)')
  patchwork::wrap_plots(panels,ncol=1)+patchwork::plot_annotation(title=paste(placement,chronotype_plot_label(gsub(' minus ',' / ',family,fixed=TRUE)),if(!is.null(context))paste(':',context)else''),
    subtitle=paste(reference,'| Factors refer to melanopic EDI + 0.1'),
    caption=paste0('A: local clock time. B: circadian phase. Factor = 10^contrast; 1 indicates no difference. Striped band: DLMO reference.\nSimultaneous bands: all 20 registered contrasts and times at 95% within each time coordinate.\nRed points exclude 1 with at least ',cfg$support_n,' participants per contributing group/day cell. Counts are unique within each signed arm.'))
}
photoperiod_figure <- function(bundle,display_data=NULL) {
  stem<-paste(bundle$label,bundle$axis,sep='_')
  if(is.null(display_data)) {
  ref<-reference_days(bundle$data);probs<-c(.1,.5,.9)
  # Weighted reference quantiles are selected from the empirical weighted CDF.
  r<-ref |> arrange(photo) |> mutate(cw=cumsum(w))
  photo_values<-vapply(probs,function(p)r$photo[which(r$cw>=p)[1]],numeric(1))
  x<-if(bundle$axis=='clock')seq(0,24-cfg$step,by=cfg$step)else seq(-12,12-cfg$step,by=cfg$step)
  V<-covariance_info(bundle$fit);z<-list()
  nd0<-prediction_template(bundle,length(x));nd0$t<-wrap24(x);nd0$photo<-photo_values[2]
  X0<-population_matrix(bundle,nd0)
  for(i in 1:3){nd<-prediction_template(bundle,length(x));nd$t<-wrap24(x);nd$photo<-photo_values[i]
    # The contrast cancels all terms except the photoperiod smooth and its time interaction.
    X<-population_matrix(bundle,nd)-X0
    z[[i]]<-bind_cols(tibble(x=x,photo=photo_values[i],reference_photo=photo_values[2],
      photoperiod=paste0(probs[i]*100,'th percentile: ',round(photo_values[i],2),' h')),linear_summary(X,coef(bundle$fit),V$V))
  }
  z<-contrast_display(bind_rows(z))
  save_table(bind_cols(z,count_sample(bundle$data)[rep(1,nrow(z)),]),paste0(stem,'_photoperiod'))
  } else z<-display_data
  p<-ggplot(z,aes(x,ratio,color=photoperiod,fill=photoperiod))
  if(bundle$axis=='phase')p<-p+dlmo_reference_band()
  p<-p+geom_hline(yintercept=1,color='grey60')+
    geom_ribbon(aes(ymin=ratio_lower,ymax=ratio_upper),alpha=.1,color=NA)+geom_line()+time_scale(bundle$axis)+factor_scale()+
    labs(title=paste(if(bundle$label=='primary')'Near-eye'else'Chest',if(bundle$axis=='phase')'circadian phase:'else'clock time:','photoperiod effects relative to the median'),subtitle=sample_label(bundle$data),x=axis_title(bundle$axis),y='Factor for (melEDI + 0.1)',
      color='Within-site centered photoperiod',fill='Within-site centered photoperiod',
      caption=paste0('Reference: the weighted median photoperiod at the same time. Bands: 95% pointwise intervals.\nFactors refer to melanopic EDI + 0.1 and use logarithmic spacing. Photoperiod is civil dawn to civil dusk.',if(bundle$axis=='phase')'\nStriped band: DLMO reference.'else''))
  save_plot(p,paste0(stem,'_photoperiod'),10,6)
}
matched_placement_check <- function(prepared) {
  key<-c('site','Id','local_date','clock_bin')
  g<-prepared$frames$glasses;c<-prepared$frames$chest
  m<-g |> inner_join(c |> select(all_of(key),chest_response=response),by=key,relationship='one-to-one') |>
    mutate(difference=chest_response-response)
  # Equal participant/day descriptive sensor difference, on exactly matched bins.
  z<-m |> summarise(mean_difference=mean(difference),RMSE=sqrt(mean(difference^2)),
    participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(quartile,daytype)) |>
    mutate(factor_at_mean_difference=10^mean_difference)
  save_table(z,'matched_placement_summary')
  hh<-m |> mutate(hour=floor(clock)) |> summarise(difference=mean(difference),.by=c(participant_key,participant_day_key,quartile,daytype,hour)) |>
    summarise(difference=mean(difference),days=n(),.by=c(participant_key,quartile,daytype,hour)) |>
    summarise(difference=mean(difference),participants=n_distinct(participant_key),participant_days=sum(days),.by=c(quartile,daytype,hour)) |>
    mutate(ratio=10^difference)
  save_table(hh,'matched_placement_hourly')
  p<-ggplot(hh,aes(hour+.5,ratio,color=quartile))+geom_hline(yintercept=1,color='grey60')+geom_line()+facet_wrap(~daytype)+factor_scale()+
    labs(title='Chest / near-eye exposure on matched half-hours',subtitle=sample_label(m),x='Local clock hour',y='Factor for (melEDI + 0.1)',
      caption='Reference: near-eye exposure. Factors are converted after averaging log differences over days within participants, then participants.\nLogarithmic axis; 1 indicates equal exposure. Counts and hourly support are exported with the figure data.')
  save_plot(p,'matched_placement',11,5)
  z
}
individual_curve_pdf <- function(curves,stem,axis,export_data=TRUE) {
  curves<-intensity_display(curves)
  if(export_data)save_table(curves,paste0(stem,'_individual_curves'))
  pdf(paste0('results/figures/',stem,'_individual_curves.pdf'),width=12,height=9)
  on.exit(dev.off())
  ids<-unique(curves$participant_key)
  for(page in split(ids,ceiling(seq_along(ids)/12))) {
    z<-filter(curves,participant_key %in% page) |> mutate(panel=paste0(participant_key,' | ',participant_days,' days'))
    p<-ggplot(z,aes(x,estimate_melEDI))
    if(axis=='phase')p<-p+dlmo_reference_band()
    print(p+geom_ribbon(aes(ymin=lower_melEDI,ymax=upper_melEDI),fill='#0072B2',alpha=.15)+
      geom_line(color='#0072B2')+geom_line(aes(y=group_estimate_melEDI),color='grey45',linetype=2)+facet_wrap(~panel,ncol=3)+
      time_scale(axis)+intensity_scale(limits=range(c(0,curves$lower_melEDI,curves$upper_melEDI)))+
      labs(title=paste(if(axis=='phase')'Circadian phase:'else'Clock time:','individual and corresponding group curves'),
        subtitle='One participant per panel; site/day effects excluded; pooled photoperiod and day-type weights',
        x=axis_title(axis),y='Melanopic EDI (lx; symlog)',
        caption=paste0('Blue: back-transformed participant curve with 95% pointwise band. Dashed: matching group curve.\nSymlog: linear below 1 lx, base-10 above 1 lx. Each panel gives the participant\'s contributing days.',if(axis=='phase')' Striped band: DLMO reference.'else'')))
  }
  invisible(curves)
}
run_outputs <- function(prepared) {
  summaries<-verdicts<-list()
  for(label in c('primary','chest')) {
    cc<-dd<-list()
    for(axis in c('clock','phase')) {
      bundle<-readRDS(paste0('results/models/',label,'_',axis,'.rds'))
      summaries[[paste(label,axis)]]<-model_summary(bundle)
      path<-paste0('results/data/',label,'_',axis,'_inference.rds')
      inference<-infer_curves(bundle)
      cc[[axis]]<-mutate(inference$curves,axis=axis);dd[[axis]]<-mutate(inference$contrasts,axis=axis)
      diagnostic<-diagnose_model(bundle);verdicts[[paste(label,axis)]]<-diagnostic$flags
      photoperiod_figure(bundle)
      if(label=='primary')individual_curves(bundle,inference)
      # Compact inference cache, retaining matrices for audit and sensitivities.
      inference$individual_draws<-NULL
      atomic_rds(inference,path)
      rm(bundle,inference);gc()
    }
    curves<-bind_rows(cc);contrasts<-bind_rows(dd)
    for(context in c('Overall','Work','Free')) {
      save_plot(pattern_figure(curves,context,if(label=='primary')'Near-eye'else'Chest'),paste0(label,'_patterns_',tolower(context)),12,11)
      for(family in c('Evening minus morning','Extreme minus moderate'))save_plot(difference_figure(contrasts,family,context,if(label=='primary')'Near-eye'else'Chest'),
        paste0(label,'_',gsub(' ','_',tolower(family)),'_',tolower(context)),12,11)
    }
    for(family in c('Free minus work','Interactions'))save_plot(difference_figure(contrasts,family,placement=if(label=='primary')'Near-eye'else'Chest'),
      paste0(label,'_',gsub(' ','_',tolower(family))),12,14)
  }
  save_table(bind_rows(summaries),'model_fit');save_table(bind_rows(verdicts),'diagnostic_verdicts')
  matched_placement_check(prepared)
  log_step('Completed primary and chest output figures')
}
assemble_figures <- function(prepared,labels=c('primary','chest')) {
  summaries<-verdicts<-list()
  for(label in labels) {
    curves<-contrasts<-list()
    for(axis in c('clock','phase')) {
      inf<-readRDS(paste0('results/data/',label,'_',axis,'_inference.rds'))
      curves[[axis]]<-mutate(inf$curves,axis=axis);contrasts[[axis]]<-mutate(inf$contrasts,axis=axis)
      stem<-paste(label,axis,sep='_')
      save_table(intensity_display(inf$curves),paste0(stem,'_curves'))
      save_table(contrast_display(inf$contrasts),paste0(stem,'_contrasts'))
      save_table(window_display(inf$windows),paste0(stem,'_significant_windows'))
      rm(inf);gc()
      b<-readRDS(paste0('results/models/',label,'_',axis,'.rds'));summaries[[paste(label,axis)]]<-model_summary(b)
      photoperiod_figure(b)
      rm(b);gc()
      if(label=='primary')individual_curve_pdf(read.csv(paste0('results/tables/',stem,'_individual_curves.csv')),stem,axis)
      verdicts[[paste(label,axis)]]<-read.csv(paste0('results/diagnostics/',label,'_',axis,'_verdict.csv'))
    }
    curves<-bind_rows(curves);contrasts<-bind_rows(contrasts)
    for(context in c('Overall','Work','Free')) {
      save_plot(pattern_figure(curves,context,if(label=='primary')'Near-eye'else'Chest'),paste0(label,'_patterns_',tolower(context)),12,11)
      for(family in c('Evening minus morning','Extreme minus moderate'))save_plot(difference_figure(contrasts,family,context,if(label=='primary')'Near-eye'else'Chest'),
        paste0(label,'_',gsub(' ','_',tolower(family)),'_',tolower(context)),12,11)
    }
    for(family in c('Free minus work','Interactions'))save_plot(difference_figure(contrasts,family,placement=if(label=='primary')'Near-eye'else'Chest'),
      paste0(label,'_',gsub(' ','_',tolower(family))),12,14)
  }
  if(length(labels)==2){save_table(bind_rows(summaries),'model_fit');save_table(bind_rows(verdicts),'diagnostic_verdicts')}
  matched_placement_check(prepared)
  log_step('Refreshed figures from saved fits and inference; no models refitted')
  invisible(NULL)
}
