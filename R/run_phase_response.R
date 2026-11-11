source('R/setup.R')
source('R/inference.R')
source('R/plots.R')
source('R/phase_response.R')
source('R/emens_balance.R')
inputs<-c('results/models/primary_phase.rds','results/data/primary_phase_inference.rds',
  'input-data/phase_response/hilaire_12_Fig3A_digitized_data.csv','R/phase_response.R',
  'input-data/phase_response/source_manifest.json','R/run_phase_response.R','R/plots.R',
  'results/tables/primary_phase_significant_windows.csv','R/emens_balance.R')
before<-tools::md5sum(inputs)
prc<-phase_response_read_prc()
bundle<-readRDS(inputs[1]);inference<-readRDS(inputs[2])
result<-phase_response_integrate(prc,bundle,inference)
phase_response_export(result)
phase_response_figures(result)
emens_balance_export(emens_balance(bundle,inference))
stopifnot(identical(before,tools::md5sum(inputs)))
save_table(tibble(path=inputs,md5=unname(before)),'phase_response_inputs','results/provenance')
capture.output(sessionInfo(),file='results/provenance/phase_response_session.txt')
writeLines(c('Rscript R/run_phase_response.R',paste('Completed',Sys.time()),
  'No primary model was refitted or changed. All scientific computation performed in R.'),
  'results/provenance/phase_response_command.txt')
