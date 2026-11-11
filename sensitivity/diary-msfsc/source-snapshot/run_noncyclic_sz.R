# Run from the project root. The fit and post-processing caches are separate
# from all primary analysis objects.
source('R/setup.R')
source('R/report_tables.R')
source('R/models.R')
source('R/inference.R')
source('R/diagnostics.R')
source('R/plots.R')
source('R/noncyclic_sz.R')
cfg$nthreads <- 1L
for(axis in c('clock','phase')) noncyclic_sz_fit(axis)
if(exists('noncyclic_sz_outputs')) noncyclic_sz_outputs()

# Execute the visible notebook's floor-diagnostic cell without maintaining a
# second copy of its analytical choices. Rendering executes this same cell.
noncyclic_clock <- readRDS('results/models/noncyclic_sz_clock.rds')
noncyclic_phase <- readRDS('results/models/noncyclic_sz_phase.rds')
notebook_lines <- readLines('noncyclic-sz-sensitivity.qmd')
cell_start <- which(notebook_lines=='#| label: tbl-sensitivity-floor')
stopifnot(length(cell_start)==1L)
cell_end <- which(seq_along(notebook_lines)>cell_start & notebook_lines=='```')[1L]
eval(parse(text=notebook_lines[seq.int(cell_start+1L,cell_end-1L)]))
