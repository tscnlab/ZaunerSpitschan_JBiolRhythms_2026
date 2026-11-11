source('R/setup.R')
source('R/inference.R')
source('R/emens_balance.R')
inputs<-c('results/models/primary_phase.rds','results/data/primary_phase_inference.rds',
  'R/emens_balance.R','R/run_emens_balance.R')
before<-tools::md5sum(inputs)
result<-emens_balance(readRDS(inputs[1]),readRDS(inputs[2]))
emens_balance_export(result)
stopifnot(identical(before,tools::md5sum(inputs)))
save_table(tibble(path=inputs,md5=unname(before)),'emens_balance_inputs','results/provenance')
writeLines(c('Command: Rscript --vanilla R/run_emens_balance.R',paste('Completed:',Sys.time()),
  paste('Seed:',result$seed,'Joint coefficient draws:',result$draws),
  'Existing fitted model and exposure curves; no refit.',
  'Source: Emens et al. 2009, doi:10.1080/07420520902821077.',
  capture.output(sessionInfo())),'results/provenance/emens_balance_runtime.txt')
print(filter(result$summary,context=='Overall',dlmo_offset== -6.45),width=Inf)
print(filter(result$contrasts,context=='Overall',contrast=='Evening / morning'),width=Inf)
