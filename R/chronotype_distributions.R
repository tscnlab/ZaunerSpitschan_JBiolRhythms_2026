# Execute the report's single source of histogram code and verify its outputs.
# Run: Rscript --vanilla R/chronotype_distributions.R
source('R/setup.R')
prepared <- readRDS('results/data/prepared.rds')
script <- tempfile('chronotype-distributions-',fileext='.R')
invisible(knitr::purl('notes/chronotype-distributions.qmd',output=script,
            documentation=0,quiet=TRUE))
sys.source(script,envir=globalenv())
unlink(script)

checks <- list()
check <- function(name,passed) {
  stopifnot(isTRUE(passed))
  checks[[length(checks)+1L]] <<- tibble(check=name,passed=TRUE)
}
check('Every near-eye participant occurs exactly once',
  nrow(distribution_people)==n_distinct(distribution_source$participant_key) &&
    !anyDuplicated(distribution_people$participant_key))
check('All participant-days and half-hours are represented',
  sum(distribution_people$participant_days)==n_distinct(distribution_source$participant_day_key) &&
    sum(distribution_people$half_hours)==nrow(distribution_source))
check('Sample is 140 participants / 778 days / 35949 half-hours',
  identical(as.integer(unlist(count_sample(distribution_source))),c(140L,778L,35949L)))
check('Neither histogram drops a missing chronotype value',
  !anyNA(distribution_people[c('msf_hours','quartile','meq','meq_type')]))
expected_quartile <- cut(distribution_people$msf_hours,
  c(-Inf,prepared$cutoffs,Inf),labels=paste0('Q',1:4),right=TRUE)
check('Frozen type-7 quartiles are retained and ties stay together',
  identical(as.character(expected_quartile),as.character(distribution_people$quartile)) &&
    isTRUE(all.equal(unname(quantile(distribution_people$msf_hours,c(.25,.5,.75),type=7)),prepared$cutoffs)))
check('Intermediate MEQ participants remain in the histogram',
  sum(distribution_people$meq_type=='Intermediate')==72L &&
    sum(distribution_people$participant_days[distribution_people$meq_type=='Intermediate'])==399L)

# Export the actual plotted bin counts and check each against individual records.
for (instrument in c('msfsc','meq')) {
  plot <- if(instrument=='msfsc')msfsc_histogram else meq_histogram
  colors <- msfsc_colors
  variable <- if(instrument=='msfsc')'msf_hours'else'meq'
  grouping <- 'quartile'
  bins <- as_tibble(ggplot_build(plot)$data[[1]]) |>
    transmute(group=names(colors)[match(fill,unname(colors))],
              xmin=xmin,xmax=xmax,midpoint=x,participants=count)
  bins <- filter(bins,is.finite(xmin),is.finite(xmax))
  # Round comparison coordinates only, below timestamp precision, so binary
  # representation of a 20-minute edge cannot exclude a point on that edge.
  comparison_values <- round(distribution_people[[variable]], 12)
  expected <- vapply(seq_len(nrow(bins)),function(i) {
    sum(as.character(distribution_people[[grouping]])==bins$group[i] &
      comparison_values>=round(bins$xmin[i], 12) &
      comparison_values<round(bins$xmax[i], 12))
  },integer(1))
  check(paste(instrument,'histogram bins agree with participant records'),
        all(bins$participants==expected))
  check(paste(instrument,'histogram contains 140 participant observations'),
        sum(bins$participants)==140L)
  check(paste(instrument,'histogram has the requested bin width'),
        all(abs((bins$xmax-bins$xmin)-if(instrument=='msfsc')20/60 else 1)<1e-12))
  save_table(bins,paste0('chronotype_distribution_',instrument,'_bins'))
}
check('Both histograms use identical MSFsc quartile colors',
      identical(ggplot_build(msfsc_histogram)$plot$scales$get_scales('fill')$palette(4),
                ggplot_build(meq_histogram)$plot$scales$get_scales('fill')$palette(4)))
check('MEQ top labels retain all five ordered category counts',
      identical(as.character(meq_category_labels$meq_type),distribution_meq_order) &&
        sum(meq_category_labels$participants)==140L &&
        sum(meq_category_labels$participant_days)==778L)
save_table(bind_rows(checks),'chronotype_distribution_acceptance','results/provenance')
files <- c('results/data/prepared.rds','notes/chronotype-distributions.qmd',
           'R/chronotype_distributions.R','R/setup.R')
save_table(tibble(file=files,md5=unname(tools::md5sum(files))),
           'chronotype_distribution_inputs','results/provenance')
writeLines(c('Command: Rscript --vanilla R/chronotype_distributions.R',
  paste('Working directory:',getwd()),paste('Executed:',Sys.time()),
  paste('OpenMP preload:',omp),paste('Passing checks:',length(checks)),
  '',capture.output(sessionInfo())),
  'results/provenance/chronotype_distribution_runtime.txt')
cat(length(checks),'chronotype distribution checks passed.\n')
