source('R/setup.R')
for(f in c('data','models','inference','diagnostics','validation','attribution','plots','sensitivity','meq','daylight')) if(file.exists(paste0('R/',f,'.R'))) source(paste0('R/',f,'.R'))
args <- commandArgs(trailingOnly=TRUE)
stage <- if(length(args))args[[1]] else 'all'
if(stage %in% c('prepare','all') || !file.exists('results/data/prepared.rds')) prepared <- prepare_data() else prepared <- readRDS('results/data/prepared.rds')
if(stage %in% c('primary','all')) {
  prep <- fit_preprocess(prepared$frames$glasses,prepared$cutoffs)
  for(axis in c('clock','phase')) { b <- fit_pair(prepared$frames$glasses,prep,axis,'primary'); print(model_summary(b)); rm(b); gc() }
}
if(stage %in% c('chest','all')) {
  prep <- fit_preprocess(prepared$frames$chest,prepared$cutoffs)
  for(axis in c('clock','phase')) { b <- fit_pair(prepared$frames$chest,prep,axis,'chest'); print(model_summary(b)); rm(b); gc() }
}
if(stage %in% c('outputs','all')) run_outputs(prepared)
if(stage=='figures') assemble_figures(prepared)
if(stage %in% c('validation','all')) run_validation(prepared)
if(stage %in% c('attribution','all')) run_attribution(prepared)
if(stage %in% c('sensitivity','all')) run_sensitivity(prepared)
if(stage %in% c('meq','all')) run_meq(prepared)
if(stage %in% c('daylight','all')) run_daylight(prepared)
if(stage %in% c('checks','all')) source('R/acceptance.R')
log_step('Completed stage',stage)
