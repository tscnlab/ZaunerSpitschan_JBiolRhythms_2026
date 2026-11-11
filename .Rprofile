# No automatic restoration or renv activation during startup.
melidos_library_root <- Sys.getenv('MELIDOS_MAIN', '../ZaunerEtAl_bioRxiv_2026')
melidos_r_series <- paste0('R-',R.version$major,'.',strsplit(R.version$minor,'.',fixed=TRUE)[[1]][1])
.libPaths(c('renv/library/standalone', Sys.getenv('MELIDOS_R_LIBRARY', ''),
  Sys.glob(file.path(melidos_library_root,'renv/library','*',melidos_r_series,R.version$platform)), .libPaths()))
melidos_omp <- Sys.getenv('MELIDOS_OPENMP_LIBRARY', '/opt/homebrew/opt/libomp/lib/libomp.dylib')
if(file.exists(melidos_omp) && !any(grepl('libomp',names(getLoadedDLLs())))) dyn.load(melidos_omp,local=FALSE)
rm(melidos_library_root,melidos_r_series,melidos_omp)
