# Run a fixed sequence in the current scenario directory, stopping on failure.
args <- commandArgs(trailingOnly = TRUE)
stages <- if (length(args)) args else
  c(
    "validation",
    "daylight",
    "daylight-validation",
    "meq",
    "sensitivity",
    "noncyclic",
    "outputs"
  )
for (stage in stages) {
  logfile <- paste0("results/logs/", stage, "_clean.log")
  # source() parses the complete driver before execution, so later edits to the
  # authoring file cannot change the running process's parsed expressions.
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "-e", shQuote('source("../scenario.R")'), shQuote(stage)),
    stdout = logfile,
    stderr = logfile
  )
  if (status != 0L) stop("Stage failed: ", stage, "; see ", logfile)
  cat(format(Sys.time()), "Completed", stage, "\n")
}
