# Run explicitly to restore the standalone analysis environment.
stopifnot(as.character(getRversion()) == '4.6.1')
if(!requireNamespace('renv',quietly=TRUE)) {
  stop('Install renv 1.2.3, then rerun Rscript --vanilla restore.R.',call.=FALSE)
}
analysis_library <- file.path(getwd(),'renv/library/standalone')
dir.create(analysis_library,recursive=TRUE,showWarnings=FALSE)
renv::restore(project=getwd(),lockfile='renv.lock',library=analysis_library,prompt=FALSE)
