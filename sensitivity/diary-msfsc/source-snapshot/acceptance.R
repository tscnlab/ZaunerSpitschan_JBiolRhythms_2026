checks<-list()
check<-function(name,condition,detail='') {
  checks[[length(checks)+1L]]<<-tibble(check=name,passed=isTRUE(condition),detail=detail)
  if(!isTRUE(condition))stop('Acceptance check failed: ',name,call.=FALSE)
}
prepared<-readRDS('results/data/prepared.rds');d<-prepared$frames$glasses
check('Primary sample',identical(as.integer(unlist(count_sample(d))),c(140L,778L,35949L)))
check('Input checksums unchanged',identical(unname(tools::md5sum(prepared$manifest$path)),prepared$manifest$md5))
check('Unique registry bins',!anyDuplicated(prepared$registry[c('placement','participant_key','local_date','clock_bin')] ))
check('Ties remain together',all((d |> summarise(groups=n_distinct(quartile),.by=msf_hours))$groups==1L))
check('Phase wraps correctly',max(abs(wrap12(d$clock-d$msf_hours)-d$phase_display))<1e-12)
check('H02 response retained',max(abs(log10(d$metric_value_lx+.1)-d$response))<1e-12 && all(d$valid_medi_wall_minutes>=15))
check('Every day starts an AR sequence',all(d$AR_start[!duplicated(d$participant_day_key)]))
check('No AR sequence crosses a gap',all(d$elapsed_from_previous_seconds[!d$AR_start]==0))
check('Ambiguous links reset AR',all(d$AR_start[!d$one_to_one_elapsed_coordinate]))
for(label in c('primary','chest')) {
  raw<-prepared$frames[[if(label=='primary')'glasses'else'chest']]
  a<-readRDS(paste0('results/models/',label,'_clock.rds'));b<-readRDS(paste0('results/models/',label,'_phase.rds'))
  check(paste(label,'identical observations and AR order'),identical(a$data$participant_day_key,b$data$participant_day_key)&&identical(a$data$clock_bin,b$data$clock_bin)&&identical(a$data$response,b$data$response)&&identical(a$data$AR_start,b$data$AR_start))
  check(paste(label,'same fixed cutoffs'),identical(a$prep$cutoffs,prepared$cutoffs)&&identical(b$prep$cutoffs,prepared$cutoffs))
  rm(a,b);gc()
  for(axis in c('clock','phase')) {
    inf<-readRDS(paste0('results/data/',label,'_',axis,'_inference.rds'))
    cr<-inf$registry
    check(paste(label,axis,'95% within-coordinate coverage'),
      identical(inf$axis_coverage,cfg$curve_coverage) &&
      abs(inf$critical-unname(quantile(inf$simulated_maxima,cfg$curve_coverage,type=8)))<1e-12)
    check(paste(label,axis,'contrast balance'),all(abs(rowSums(cr$W))<1e-10)&&nrow(cr$W)==20)
    check(paste(label,axis,'evening-minus-morning sign'),cr$W['evening_morning_Extreme_Work','Q4_Work']==1&&cr$W['evening_morning_Extreme_Work','Q1_Work']==-1)
    check(paste(label,axis,'support gate'),!any(inf$contrasts$significant & !inf$contrasts$supported))
    for(severity in c('Moderate','Extreme'))for(context in c('Overall','Work','Free')) {
      qs<-if(severity=='Moderate')c('Q2','Q3')else c('Q1','Q4')
      lhs<-filter(inf$curves,context==.env$context,quartile==qs[2]) |> arrange(x)
      rhs<-filter(inf$curves,context==.env$context,quartile==qs[1]) |> arrange(x)
      contrast_id<-paste('evening_morning',severity,context,sep='_')
      difference<-filter(inf$contrasts,contrast==.env$contrast_id) |> arrange(x)
      check(paste(label,axis,severity,context,'contrast equals plotted curve difference'),nrow(lhs)==nrow(difference)&&nrow(rhs)==nrow(difference)&&max(abs(lhs$estimate-rhs$estimate-difference$estimate))<1e-9)
    }
    for(context in c('Overall','Work','Free'))for(q in paste0('Q',1:4)) {
      rows<-filter(raw,as.character(quartile)==q)
      if(context!='Overall')rows<-filter(rows,as.character(daytype)==context)
      count<-count_sample(rows);shown<-filter(inf$curves,.data$context==.env$context,quartile==q)
      check(paste(label,axis,q,context,'displayed sample'),all(shown$participants==count$participants)&all(shown$participant_days==count$participant_days)&all(shown$half_hours==count$half_hours))
    }
    for(j in seq_len(nrow(inf$windows))) {
      win<-inf$windows[j,];w<-cr$W[win$contrast,];coord<-if(axis=='clock')raw$clock else raw$phase_display
      rows<-raw[wrap24(coord-win$start)<win$duration_hours-1e-10,]
      pos<-counts_for_cells(rows,names(w)[w>0]);neg<-counts_for_cells(rows,names(w)[w<0])
      check(paste(label,axis,'window',j,'sample'),pos$participants==win$positive_participants&&pos$participant_days==win$positive_days&&pos$half_hours==win$positive_half_hours&&neg$participants==win$negative_participants&&neg$participant_days==win$negative_days&&neg$half_hours==win$negative_half_hours)
    }
    rm(inf);gc()
  }
}
for(axis in c('clock','phase')) {
  sh<-read.csv(paste0('results/tables/shapley_',axis,'.csv'))
  check(paste(axis,'Shapley sum'),abs(sum(sh$allocation)-(sh$full_R2[1]-sh$baseline_R2[1]))<1e-9)
  subsets<-read.csv(paste0('results/tables/shapley_subsets_',axis,'.csv'))
  check(paste(axis,'all subsets'),identical(subsets$mask,0:127)&&all(subsets$rho==subsets$rho[1]))
}
splits<-validation_splits(d)
for(s in splits){train<-d[!s$assessment,];test<-d[s$assessment,]
  check(paste(s$id,'no day overlap'),length(intersect(train$participant_day_key,test$participant_day_key))==0)
  if(s$target=='New participants')check(paste(s$id,'no participant overlap'),length(intersect(train$participant_key,test$participant_key))==0)
  for(axis in c('clock','phase')){
    v<-readRDS(paste0('results/validation/',s$id,'_',axis,'.rds'))
    check(paste(s$id,axis,'training cutoffs'),identical(v$cutoffs,quartile_cutoffs(train)))
    check(paste(s$id,axis,'assessment counts'),identical(as.integer(unlist(v$assessment)),as.integer(unlist(count_sample(test)))))
  }
}
pred<-read.csv('results/validation/held_out_predictions.csv')
ca<-filter(pred,axis=='clock') |> arrange(target,participant_key,local_date,clock_bin)
pa<-filter(pred,axis=='phase') |> arrange(target,participant_key,local_date,clock_bin)
check('Matched predictive assessment rows',identical(ca$observed,pa$observed)&&identical(ca$participant_day_key,pa$participant_day_key)&&identical(ca$clock_bin,pa$clock_bin))
check('DLMO guide endpoints',max(abs(dlmo_reference()-c(-6.95,-5.95)))<1e-10)
check('Circular boundary merge represented',any(read.csv('results/tables/primary_clock_significant_windows.csv')$wraps_boundary))
save_table(bind_rows(checks),'acceptance_checks','results/provenance')
sources<-c(list.files('R',full.names=TRUE,pattern='\\.R$'),list.files(pattern='\\.(qmd|yml|bib|R|css)$'),'notes/chronotype-distributions.qmd')
save_table(tibble(path=sources,md5=unname(tools::md5sum(sources))),'analysis_sources','results/provenance')
log_step('All',length(checks),'acceptance checks passed')
