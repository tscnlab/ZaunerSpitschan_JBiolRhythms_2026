# The everyday light exposome

This Quarto project reproduces the manuscript **The health-relevant architecture of the everyday light exposome**, its supplementary information, and the accompanying preparation and analysis pages. The analysis uses the local data supplied in this repository.

> **Current and legacy analyses.** The `main` branch and [current analysis website](https://tscnlab.github.io/ZaunerEtAl_bioRxiv_2026/) contain the latest analysis, prepared for [bioRxiv version 3](https://www.biorxiv.org/content/10.64898/2026.06.04.730226v3), whose update is forthcoming. The analysis reported in [version 1](https://www.biorxiv.org/content/10.64898/2026.06.04.730226v1) and [version 2](https://www.biorxiv.org/content/10.64898/2026.06.04.730226v2) is preserved in the [legacy branch](https://github.com/tscnlab/ZaunerEtAl_bioRxiv_2026/tree/legacy).

## Software and environment

Use R 4.6.1 and Quarto 1.9.37 for the reference environment. R package versions are recorded in `renv.lock`. A working C++ compiler for R is required because the recommendation models compile their likelihood and prediction code with TMB and Rcpp. On macOS this requires the Xcode Command Line Tools; on Windows use the Rtools release appropriate for your R version.

Open `reproducible-NH.Rproj` in RStudio, or open a terminal in the repository root. Restore the recorded R packages once:

```sh
Rscript -e 'renv::restore(prompt = FALSE); renv::status()'
```

Package restoration needs internet access unless the packages are already cached. The analysis reads the supplied local data and does not require downloading the study datasets again. Do not update packages or take a new snapshot before comparing a reproduction with this baseline.

The project uses `config/Makevars` when compiling R packages and model code, so unrelated personal compiler settings do not change the build. It uses R's compiler defaults and package-specific configuration. On macOS it also discovers an existing Homebrew gettext installation, whose development files may be needed to build `nlme` from source. Package cache reuse is disabled for this project to avoid reusing binaries built with incompatible compiler settings; package versions remain fixed by `renv.lock`.

The project disables automatic dependency scanning during R startup. An explicit `renv::status()` still checks the environment; `.renvignore` excludes downloaded data, generated results, and rendered pages from dependency discovery. RStudio starts without restoring or saving an `.RData` workspace.

## Run the complete project

From the repository root, run:

```sh
quarto render
```

This command executes preparation and analyses in the order declared in `_quarto.yml`, then assembles the manuscript and supplementary information. It regenerates figures, tables, numerical summaries, fitted models, the website in `docs/`, and the manuscript Word file. No earlier fitted models or rendered pages are required. Each analysis page explains its inputs, analytical choices, and outputs alongside executable R code. The website menus provide separate links to the data-and-model guide and the findings within each analysis page.

For a shorter run with **50 bootstrap replicates** wherever bootstrap sampling is used:

```sh
quarto render --profile quick
```

The full run retains each analysis's stated default count. The quick profile changes bootstrap precision, so its intervals are for checking execution and should not replace the full manuscript results. Model fitting and data preparation still run in full. Simulation-based diagnostic checks also use the central quick count. The central quick count and available parallel workers are set in `config/analysis.yml`. Table operations use one data.table thread; model and bootstrap workers are controlled separately. Both commands overwrite the same generated output locations; finish with a full run before interpreting or distributing results.

On the reference machine, an Apple M4 Max with 16 CPU cores and 128 GiB of memory, using four workers, R 4.6.1, and Quarto 1.9.37, the complete render from source inputs took approximately **8 hours 54 minutes**. Package restoration was completed before timing the render.

| Run | Approximate wall-clock time on the reference machine |
|---|---|
| `quarto render` | **8 hours 54 minutes**, measured |
| `quarto render --profile quick` (50 bootstrap replicates) | **5–6 hours**, estimated |

The quick-run estimate is based on the measured full run and a timed reduced-bootstrap run of the site-differences analysis; it is not an end-to-end quick-run benchmark. Reducing bootstrap counts saves substantial time in that analysis, but preparation, model fitting, and many sensitivity fits still execute in full. Allow for variation with hardware, available workers, and other activity on the machine.

## Repository contents

```text
data/
  downloaded/                 Fixed recordings, questionnaires, and diaries
  alternative-baseline/       Local inputs for alternative preprocessing
preparation/                  Import, alignment, coverage, metrics, and samples
analyses/                     Models, sensitivities, summaries, and displays
scripts/                      Shared preparation, model, and display functions
config/                       Data-source definitions and analysis settings
assets/                       Study diagram and Word formatting template
styles/                       Manuscript styling
_extensions/                  Required Quarto formatting code
index.qmd                     Manuscript
supplementary.qmd             Supplementary information entry point
_supplementary.qmd            Supplementary content included by that page
preregistration-deviations.qmd Scientific deviations and their rationale
references.bib                Bibliography
nature.csl                    Citation style
_quarto.yml                   Project settings and execution order
_quarto-quick.yml             Reduced-bootstrap profile
renv.lock, renv/, .Rprofile    Reproducible R environment
reproducible-NH.Rproj          RStudio project
```

The alternative-preprocessing analysis starts from the six supplied `.RData` files in `data/alternative-baseline/`. These are fixed data inputs, not fitted model outputs. The preparation and alternative-preprocessing pages explain which temporal and support rules differ and how the sensitivity datasets are derived.

The shared scripts define functions used by the Quarto pages. They do not need separate terminal commands. Analysis settings and the calls that perform each calculation remain in the relevant page. The Word layout hook is the only post-render step: `scripts/reporting/word-layout.R` keeps tables editable, repeats their headers, prevents split rows, and places wide tables on landscape pages. Quarto runs it automatically; no separate command or manual table correction is required. The render creates the manuscript at `docs/index.docx` and a separate supplementary methods document at `docs/preregistration-deviations.docx`. The hook also inserts section and supplementary-display page breaks and normalizes Word table markup. Close the manuscript in Word before rendering it again. If Word reports an editing lock after an interrupted session even though no other copy is open, quit Word normally and reopen it; the generated document has no editing protection.

## Generated results

Rendering creates the following output directories:

```text
docs/                         Rendered website and manuscript Word document
results/
  intermediate/               Prepared measurements and analysis datasets
  models/                     Fitted model objects and fitted samples
  images/                     Exported figures
  tables/                     Editable tables and their numerical data
  csv/
    source_data/              Figure data, estimates, and summaries
    diagnostics/              Numerical model and support checks
```

The tables and CSV files support checking numerical results without reading values from plots. Random sampling uses explicit seeds. Small numerical differences may arise from different operating systems, numerical libraries, or compiler versions; compare the stated estimands and uncertainty as well as individual coefficients.

Git tracks the rendered website in `docs/`, preprocessing outputs in `results/intermediate/` (including RDS files), and the exported figures, tables, and CSV results. Serialized analysis objects, identified by RDS files elsewhere under `results/`, are omitted from Git for size reasons: saved model objects can exceed [GitHub's 100 MiB per-file limit](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github) and substantially enlarge the repository. These files remain available locally and are regenerated by the analysis pages during `quarto render`; they are not required to serve the rendered website. Generated model compilation files are also excluded and rebuilt locally. The `.renvignore` file separately excludes all generated outputs from dependency discovery, whether or not Git tracks them.

For GitHub Pages, publish the committed `docs/` directory. Quarto copies the root `.nojekyll` marker into that directory so the rendered site can be served without additional Jekyll processing. After a full render, review and commit updated source files together with the generated website and tracked results.

To reproduce from a source-only starting point, use a separate checkout and remove `docs/`, `results/`, `.quarto/`, and any root `site_libs/` directory there. Quarto may also leave generated `.html` or `.docx` files beside the pages listed in `_quarto.yml`, together with `_files/`, `_cache/`, `.rmarkdown`, or `.knit.md` intermediates; these can be removed too. Keep `data/`, `renv/`, `renv.lock`, the source files, and static assets such as `assets/reference.docx`. A subsequent `quarto render` reconstructs the outputs. Removing tracked outputs appears as deletions in Git until they are regenerated.
