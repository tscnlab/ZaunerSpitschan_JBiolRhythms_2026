# Exploratory projection of the existing phase curves onto a published PRC.
# Sourcing this file defines functions only. The report calls each calculation.
phase_response_read_prc <- function(path='input-data/phase_response/hilaire_12_Fig3A_digitized_data.csv') {
  stopifnot(file.exists(path))
  z<-read.csv(path)
  stopifnot(all(c('initial_phase','phase_shift','dlmo_collection_method') %in% names(z)),
            nrow(z)==28L,all(is.finite(z$initial_phase)),all(is.finite(z$phase_shift)))
  z$source_row<-seq_len(nrow(z))
  # This is a previously digitized subject-level scatter plot, deposited by
  # Walch et al. It is not original raw data. Fig. 3A is double plotted.
  used<-z[z$initial_phase>= -12 & z$initial_phase<12,]
  extra<-z[z$initial_phase< -12 | z$initial_phase>=12,]
  stopifnot(nrow(used)==18L,nrow(extra)==10L)
  duplicate_audit<-bind_rows(lapply(seq_len(nrow(extra)),function(i){
    delta<-abs(wrap12(used$initial_phase-extra$initial_phase[i]));j<-which.min(delta)
    tibble(excluded_source_row=extra$source_row[i],retained_source_row=used$source_row[j],
      circular_time_difference=delta[j],shift_difference=abs(used$phase_shift[j]-extra$phase_shift[i]),
      same_collection_method=used$dlmo_collection_method[j]==extra$dlmo_collection_method[i])
  }))
  stopifnot(all(duplicate_audit$circular_time_difference<.05),
    all(duplicate_audit$shift_difference<.01),all(duplicate_audit$same_collection_method))
  # A 1 h pulse onset becomes its midpoint by adding 0.5 h. Fit the same
  # two-harmonic form as St Hilaire et al. to the 18 retained coordinates.
  used$hours_from_dlmo_midpoint<-used$initial_phase+.5
  used$phase_shift_hours<-used$phase_shift
  basis<-function(x)cbind(Intercept=1,cos1=cos(2*pi*x/24),sin1=sin(2*pi*x/24),
    cos2=cos(4*pi*x/24),sin2=sin(4*pi*x/24))
  fit<-lm.fit(basis(used$hours_from_dlmo_midpoint),used$phase_shift_hours)
  rss<-sum(fit$residuals^2);tss<-sum((used$phase_shift_hours-mean(used$phase_shift_hours))^2)
  extrema<-range(drop(basis(seq(-12,12,length.out=24001))%*%fit$coefficients))
  # External checks use the paper's rounded fitted extrema, not a residual
  # threshold: residual scatter represents actual between-person variation.
  stopifnot(abs(extrema[1]-(-1.75))<.02,abs(extrema[2]-.45)<.02)
  list(points=used,all_points=z,duplicate_audit=duplicate_audit,
    coefficients=fit$coefficients,basis=basis,
    predict=function(x)drop(basis(x)%*%fit$coefficients),
    fit_rmse=sqrt(mean(fit$residuals^2)),r_squared=1-rss/tss,
    adjusted_r_squared=1-(rss/(nrow(used)-5))/(tss/(nrow(used)-1)),
    minimum_shift=extrema[1],maximum_shift=extrema[2],source_file=path)
}

phase_response_weights <- function(prc,x,dlmo_offset=-6.45,correct_drift=TRUE) {
  # x is a bin MIDPOINT relative to MSFsc. The PRC abscissa is the light
  # pulse MIDPOINT relative to DLMO (converted from onset on import).
  kernel<-prc$predict(x-dlmo_offset)+if(correct_drift).54 else 0
  step<-24/length(x)
  denominator<-sum(abs(kernel))*step
  stopifnot(denominator>0)
  tibble(x=x,dlmo_offset=dlmo_offset,
    prc_baseline=if(correct_drift)'Drift corrected'else'Published raw',
    kernel=kernel,weight=kernel*step/denominator,step=step,
    normalization=denominator)
}

phase_response_integrate <- function(prc,bundle,inference) {
  stopifnot(bundle$axis=='phase',identical(names(inference$pattern_matrices),
    paste(rep(paste0('Q',1:4),3),rep(c('Overall','Work','Free'),each=4),sep='_')))
  x<-inference$cell_design$x
  X<-inference$pattern_matrices[grepl('^Q[14]_',names(inference$pattern_matrices))]
  stopifnot(length(X)==6L,length(x)==144,all(abs(x-seq(-12,12-1/6,by=1/6))<1e-10))
  beta<-coef(bundle$fit);cv<-covariance_info(bundle$fit)
  keys<-names(X);keep<-which(Reduce(`+`,lapply(X,function(m)colSums(abs(m))))>1e-12)
  XX<-lapply(X,function(m)m[,keep,drop=FALSE]);b<-beta[keep];V<-cv$V[keep,keep,drop=FALSE]
  curves<-lapply(XX,function(m)drop(m%*%b))
  counts<-inference$curves |> filter(quartile %in% c('Q1','Q4')) |>
    distinct(quartile,context,direction,extremity,participants,participant_days,half_hours) |>
    mutate(key=paste(quartile,context,sep='_'))
  scales<-'Zero-exposure centered'
  alignments<-c(-8.45,-6.95,-6.45,-5.95,-4.45)
  results<-components<-contrasts<-weights<-checks<-list()
  add_check<-function(name,passed,value=NA_real_)checks[[length(checks)+1L]]<<-
    tibble(check=name,passed=isTRUE(passed),value=value)
  for(drift in c(TRUE,FALSE))for(offset in alignments) {
    w<-phase_response_weights(prc,x,offset,drift);weights[[length(weights)+1L]]<-w
    signed_mass<-sum(w$weight);integrated<-lapply(XX,function(m)drop(crossprod(w$weight,m)))
    add_check(paste('Absolute weights sum to one',drift,offset),abs(sum(abs(w$weight))-1)<1e-12)
    for(key in keys) {
      a<-integrated[[key]];raw<-sum(a*b);se<-sqrt(max(0,drop(a%*%V%*%a)))
      eta<-curves[[key]];row<-counts[counts$key==key,]
      for(base in scales) {
        reference<- -1
        point<-raw-reference*signed_mass;ci<-point+c(-1,1)*qnorm(.975)*se
        results[[length(results)+1L]]<-bind_cols(row,tibble(dlmo_offset=offset,
          reference_alignment=offset>=-6.95 & offset<=-5.95,
          prc_baseline=w$prc_baseline[1],exposure_scale=base,reference_log10=reference,
          reference_melEDI=10^reference-.1,estimate=point,se=se,lower=ci[1],upper=ci[2],
          unnormalized_integral=point*w$normalization[1],
          unnormalized_lower=ci[1]*w$normalization[1],unnormalized_upper=ci[2]*w$normalization[1],
          normalization=w$normalization[1],kernel_signed_mass=signed_mass,
          index_direction=if(point>0)'Advance-weighted'else'Delay-weighted',
          interval_direction=if(ci[1]>0)'Advance-weighted'else if(ci[2]<0)'Delay-weighted'else'Includes zero',
          interval_method='Integrated design; Gaussian coefficient interval; zero-exposure reference',
          covariance=cv$label))
        error<-abs(point-sum((eta-reference)*w$weight))
        add_check(paste('Direct integral agrees',drift,offset,key,base),error<1e-10,error)
        if(offset== -6.45) components[[length(components)+1L]]<-bind_cols(
          tibble(key=key,x=x,prc_baseline=w$prc_baseline[1],exposure_scale=base,
            reference_log10=reference,fitted_log10=eta,centered_log10=eta-reference,
            kernel=w$kernel,weight=w$weight,contribution=(eta-reference)*w$weight),
          row[rep(1,length(x)),setdiff(names(row),'key')])
      }
    }
    for(context in c('Overall','Work','Free')) {
      kk<-paste(c('Q1','Q4'),context,sep='_');a<-integrated[[kk[2]]]-integrated[[kk[1]]]
      raw<-sum(a*b);se<-sqrt(max(0,drop(a%*%V%*%a)))
      a_count<-counts[counts$key==kk[2],];b_count<-counts[counts$key==kk[1],]
      for(base in scales) {
        ref_difference<-0
        point<-raw-ref_difference*signed_mass;ci<-point+c(-1,1)*qnorm(.975)*se
        contrasts[[length(contrasts)+1L]]<-tibble(context=context,extremity='Extreme',
          contrast='Evening minus Morning',dlmo_offset=offset,
          reference_alignment=offset>=-6.95 & offset<=-5.95,
          prc_baseline=w$prc_baseline[1],exposure_scale=base,
          reference_log10_difference=ref_difference,estimate=point,se=se,lower=ci[1],upper=ci[2],
          evening_participants=a_count$participants,evening_days=a_count$participant_days,
          evening_half_hours=a_count$half_hours,morning_participants=b_count$participants,
          morning_days=b_count$participant_days,morning_half_hours=b_count$half_hours)
      }
    }
  }
  add_check('24 h phase grid',abs(length(x)/6-24)<1e-12)
  add_check('Published PRC cyclic continuity',abs(prc$predict(-12)-prc$predict(12))<1e-10)
  add_check('Source retained one period and 18 participants',nrow(prc$points)==18L,nrow(prc$points))
  add_check('Source duplicates match',all(prc$duplicate_audit$circular_time_difference<.05 &
    prc$duplicate_audit$shift_difference<.01 & prc$duplicate_audit$same_collection_method))
  add_check('Fitted minimum agrees with published rounded value',abs(prc$minimum_shift+1.75)<.02,prc$minimum_shift)
  add_check('Fitted maximum agrees with published rounded value',abs(prc$maximum_shift-.45)<.02,prc$maximum_shift)
  z<-list(summary=bind_rows(results),contrasts=bind_rows(contrasts),components=bind_rows(components),
    weights=bind_rows(weights),
    checks=bind_rows(checks),covariance=cv$label,prc=prc,sample=count_sample(bundle$data),
    extreme_sample=count_sample(filter(bundle$data,quartile %in% c('Q1','Q4'))))
  stopifnot(nrow(z$summary)==60L,nrow(z$contrasts)==30L,all(z$checks$passed))
  z
}

phase_response_export <- function(result) {
  for(n in c('summary','contrasts','components','weights','checks'))
    save_table(result[[n]],paste0('phase_response_',n),
      dir=if(n=='checks')'results/provenance'else'results/tables')
  save_table(tibble(term=names(result$prc$coefficients),coefficient=unname(result$prc$coefficients)),
    'phase_response_fitted_coefficients','results/provenance')
  save_table(result$prc$points,'phase_response_retained_source_points','results/provenance')
  save_table(result$prc$duplicate_audit,'phase_response_duplicate_audit','results/provenance')
  save_table(tibble(fit_rmse=result$prc$fit_rmse,r_squared=result$prc$r_squared,
    adjusted_r_squared=result$prc$adjusted_r_squared,minimum_shift=result$prc$minimum_shift,
    maximum_shift=result$prc$maximum_shift,source_participants=nrow(result$prc$points),
    source=result$prc$source_file),'phase_response_source_fit','results/provenance')
  save_table(result$summary |> filter(context=='Overall',dlmo_offset== -6.45) |>
      select(quartile,direction,extremity,participants,participant_days,half_hours,
        prc_baseline,exposure_scale,reference_log10,reference_melEDI,estimate,lower,upper,unnormalized_integral),
    'phase_response_overall_report')
  save_table(result$contrasts |> filter(dlmo_offset== -6.45),'phase_response_contrast_report')
  atomic_rds(result,'results/data/phase_response_results.rds')
  invisible(result)
}

phase_response_kernel_figure <- function(prc) {
  offset<- -6.45;grid<-seq(-12,12,by=1/60)
  z<-tibble(phase_msfsc=grid,published=prc$predict(grid-offset),drift_corrected=published+.54) |>
    pivot_longer(c(published,drift_corrected),names_to='curve',values_to='response')
  points<-as_tibble(prc$points) |> mutate(phase_msfsc=wrap12(hours_from_dlmo_midpoint+offset))
  labels<-c(published='Published raw: P(u)',drift_corrected='Drift corrected: P(u) + 0.54 h')
  ylimits<-range(c(z$response,points$phase_shift_hours))+c(-.12,.65)
  p<-ggplot(z,aes(phase_msfsc,response,colour=curve,linetype=curve))+
    primary_phase_interval_band()+dlmo_reference_band()+geom_hline(yintercept=0,colour='grey55',linewidth=.45)+
    geom_vline(xintercept=offset,colour='#8B6423',linewidth=.55,linetype=3)+
    geom_vline(xintercept=0,colour='#0072B2',linewidth=.6,linetype=3)+
    geom_point(data=points,aes(phase_msfsc,phase_shift_hours),inherit.aes=FALSE,colour='grey50',alpha=.75,size=1.7)+
    geom_line(linewidth=.85)+
    annotate('label',x=offset,y=ylimits[2]-.12,label='DLMO reference\n-6.45 h',colour='#8B6423',fill='white',linewidth=0,vjust=1,size=3.5)+
    annotate('label',x=0,y=ylimits[2]-.12,label='MSFsc\n0 h',colour='#0072B2',fill='white',linewidth=0,vjust=1,size=3.5)+
    scale_colour_manual(values=c(drift_corrected='#6A3D9A',published='#666666'),breaks=names(labels),labels=labels)+
    scale_linetype_manual(values=c(drift_corrected='solid',published='longdash'),breaks=names(labels),labels=labels)+
    scale_x_continuous(breaks=seq(-12,12,3),limits=c(-12,12),
      sec.axis=sec_axis(~.-offset,breaks=seq(-6,18,6),name='Light-pulse midpoint: hours relative to DLMO'))+
    scale_y_continuous(limits=ylimits,expand=expansion(mult=c(.01,.01)))+
    labs(title='Light phase-response curve anchored to DLMO and MSFsc',
      subtitle='Central alignment. Positive responses indicate advances; negative responses indicate delays.',
      x='Light-pulse midpoint: hours relative to MSFsc',y='Phase response (hours)',colour=NULL,linetype=NULL,
      caption=paste0('St Hilaire et al. (2012); numerical coordinates deposited by Walch et al. (2025). Grey points: 18 retained source coordinates.\n',
        'Both curves are shown before weight normalization. DLMO is aligned 6.45 h before MSFsc.\n',
        'Grey windows with dashed boundaries: primary MeLiDos phase intervals.
Striped band (-6.95 to -5.95 h): published group-average DLMO guide, not individual measurements.'))+
    guides(colour=guide_legend(nrow=1),linetype=guide_legend(nrow=1))
  save_plot(p,'phase_response_kernel',12.5,7)
  invisible(p)
}

phase_response_figures <- function(result) {
  z<-result$summary |> filter(dlmo_offset== -6.45) |>
    mutate(prc_baseline=factor(prc_baseline,c('Published raw','Drift corrected')),
      context=factor(context,c('Overall','Work','Free')),
      group_direction=factor(direction,c('Morning','Evening')))
  counts<-z |> distinct(context,group_direction,participants,participant_days) |>
    arrange(context,group_direction) |>
    summarise(label=paste(paste0(group_direction,': ',participants,' participants / ',participant_days,' days'),collapse='; '),.by=context)
  p<-ggplot(z,aes(context,estimate,colour=group_direction,group=group_direction))+
    geom_hline(yintercept=0,colour='grey55',linetype=2)+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=.15,position=position_dodge(.45))+
    geom_point(position=position_dodge(.45),size=2.5)+
    facet_wrap(~prc_baseline,nrow=1)+
    scale_colour_manual(values=c(Morning='#0072B2',Evening='#D55E00'))+
    labs(x=NULL,y='Signed PRC-weighted index',colour='Extreme chronotype',
      title='PRC weighting in extreme morning and evening groups',
      subtitle='Positive: advance-weighted; negative: delay-weighted. Indices are not hours of phase shift.',
      caption=paste0('Central DLMO reference: 6.45 h before MSFsc. Zero-exposure centered. 95% coefficient intervals, conditional on PRC and alignment.\n',
        paste(paste0(counts$context,': ',counts$label),collapse='\n'),'\n',
        'Full fitted model: ',result$sample$participants,' participants; ',result$sample$participant_days,' days; ',
        format(result$sample$half_hours,big.mark=','),' half-hours. Extreme groups: ',
        result$extreme_sample$participants,' participants; ',result$extreme_sample$participant_days,' days; ',
        format(result$extreme_sample$half_hours,big.mark=','),' half-hours.'))+
    theme(strip.text.y=element_text(size=10))
  save_plot(p,'phase_response_indices',12,5.8)
  phase_response_kernel_figure(result$prc)
  invisible(result)
}
