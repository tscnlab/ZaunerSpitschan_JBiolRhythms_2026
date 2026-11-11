# Verify display changes and the new Emens-style sensitivity from saved outputs.
source('R/setup.R')
checks<-list()
check<-function(name,passed,value=NA_real_)checks[[length(checks)+1L]]<<-
  tibble(check=name,passed=isTRUE(passed),value=value)
protected<-read.csv('results/provenance/revision9_protected_before.csv')
check('Existing models, PRC results and daylength results preserved',
  all(unname(tools::md5sum(protected$path))==protected$md5))
plot_curves<-read.csv('data/intermediate/manuscript/figure-2-patterns-curves.csv')
for(axis in c('clock','phase')) {
  inf<-readRDS(paste0('results/data/primary_',axis,'_inference.rds'))
  check(paste('Exposure simultaneous family',axis),inf$pattern_coverage==.95 &&
    length(inf$pattern_simulated_maxima)==10000L && all(inf$curves$band_type=='simultaneous'))
  check(paste('Exposure critical matches simulated quantile',axis),
    abs(inf$pattern_critical-quantile(inf$pattern_simulated_maxima,.95,type=8))<1e-12)
  tab<-read.csv(paste0('results/tables/primary_',axis,'_curves.csv')) |>
    filter(context=='Overall') |> arrange(quartile,x)
  p<-plot_curves |> filter(.data$axis==.env$axis) |> arrange(quartile,x)
  cols<-c('estimate_melEDI','lower_melEDI','upper_melEDI')
  check(paste('Figure 2 exposure intervals match exports',axis),
    isTRUE(all.equal(as.matrix(p[,cols]),as.matrix(tab[,cols]),tolerance=1e-10)))
}
emens<-readRDS('results/data/emens_balance.rds')
check('Emens internal region and denominator checks pass',all(emens$checks$passed))
check('Emens reports only extreme groups',setequal(emens$summary$quartile,c('Q1','Q4')))
check('No point estimate needs zero flooring',all(emens$floor_counts$negative_fitted_points==0))
for(offset in unique(emens$regions$dlmo_offset)) {
  r<-filter(emens$regions,dlmo_offset==offset)
  # Reconstruct the advance region directly as [-16,-4) relative to DLMO.
  u<-((r$x-offset+16)%%24)-16
  check(paste('Advance mask matches published limits',offset),
    all((r$region=='Advance')==(u>= -16 & u< -4)))
  check(paste('Each region spans 12 hours',offset),all(table(r$region)==72L))
}
g<-filter(emens$summary,dlmo_offset== -6.45,context=='Overall') |> arrange(quartile)
c<-filter(emens$contrasts,dlmo_offset== -6.45,context=='Overall',contrast=='Evening / morning')
check('Direct comparison equals ratio of fitted balances',abs(c$estimate-g$estimate[2]/g$estimate[1])<1e-12)
check('All group ratios at guide alignments exceed one with individual intervals',
  all(filter(emens$contrasts,contrast=='Evening / morning')$lower>1))
p<-read.csv('data/intermediate/manuscript/figure-4-emens-balances.csv') |>
  arrange(quartile,context)
expected<-filter(emens$summary,dlmo_offset== -6.45) |> arrange(quartile,context)
cols<-c('estimate','lower','upper')
check('Figure 4C exactly matches Emens estimates and intervals',
  isTRUE(all.equal(as.matrix(p[,cols]),as.matrix(expected[,cols]),tolerance=1e-10)))
z<-bind_rows(checks)
save_table(z,'revision9_verification','results/provenance')
writeLines(c('Command: Rscript --vanilla results/provenance/revision9_verify.R',
  paste('Executed:',Sys.time()),capture.output(sessionInfo())),
  'results/provenance/revision9_verification_runtime.txt')
print(z,width=Inf)
stopifnot(all(z$passed))
