<img src="https://github.com/tscnlab/Templates/blob/main/logo/logo_with_text-01.png" alt="Translational Sensory and Circadian Neuroscience" width="400"/>

# Light Exposure and Chronotype

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24140946.svg)](https://doi.org/10.5281/zenodo.24140946)

This repository accompanies the working manuscript **Daily light exposure differences between extreme chronotypes persist after alignment to habitual sleep timing**, by Johannes Zauner and Manuel Spitschan.

[Read the manuscript](https://tscnlab.github.io/ZaunerSpitschan_JBiolRhythms_2026/) · [Online Supplement S1](https://tscnlab.github.io/ZaunerSpitschan_JBiolRhythms_2026/analysis.html) · [Online Supplement S2](https://tscnlab.github.io/ZaunerSpitschan_JBiolRhythms_2026/diary-msfsc-sensitivity.html)

## Scope

The analysis compares daily melanopic light exposure across four chronotype groups in the MeLiDos cohort. The primary near-eye sample contains 140 participants, 778 participant-days and 35,949 half-hour observations. Cyclic generalized additive mixed models compare exposure in local clock time and relative to questionnaire-derived midsleep on free days corrected for sleep debt (MSFsc).

The project also examines work/free-day differences, photoperiod moderation, questionnaire and sensor-placement sensitivities, prediction and attribution, and agreement between sleep-diary timing and MSFsc. Exploratory PRC weighting is a conceptual summary of exposure timing. It is not calibrated to predict phase shifts. A separate sensitivity repeats the analyses using diary-derived corrected midsleep, with matched-sample and fixed-classification comparisons.

The manuscript and supplements are drafts and will be updated. The [Zenodo archive](https://doi.org/10.5281/zenodo.24140946) provides a persistent DOI for the repository across releases. For exact reproducibility, record the version-specific DOI or Git commit used.

## Project layout

| Path | Purpose |
|---|---|
| `index.qmd` | Manuscript, with HTML and Word formats |
| `analysis.qmd` | Online Supplement S1: primary analysis and supporting checks |
| `diary-msfsc-sensitivity.qmd` | Online Supplement S2: diary-derived timing sensitivity |
| `noncyclic-sz-sensitivity.qmd` | Additional spline-boundary sensitivity |
| `figures.qmd` | Manuscript and standalone figure downloads |
| `R/`, `run.R` | Analytical functions and stage driver |
| `sensitivity/diary-msfsc/` | Isolated sensitivity drivers and comparison exports |
| `input-data/` | Fixed prepared source inputs, provenance and upstream license |
| `results/`, `data/intermediate/` | Saved figures, tables and supporting exports |
| `_freeze/` | Quarto computational outputs used for website rendering |
| `scripts/`, `styles/` | Publishing checks, Word formatting and reference template |
| `renv.lock`, `restore.R` | Exact R package records and environment restoration |

Large fitted models, simulation caches, local package libraries, rendered site files and editorial review documents are excluded from Git. The source code and fixed prepared inputs are retained to regenerate the analyses. The website is built from committed computational outputs, so publication does not repeat model fitting.

## Render the website from saved outputs

Install **Quarto 1.9.37**, clone the repository, then run a whole-project render:

```sh
git clone https://github.com/tscnlab/ZaunerSpitschan_JBiolRhythms_2026.git
cd ZaunerSpitschan_JBiolRhythms_2026
quarto render --to html
python3 scripts/check-site.py report
quarto preview --no-browser
```

The site is written to `report/`. The committed `_freeze/` directory supplies the results of R cells. A whole-project render reuses these results when the computational documents have not changed. Rendering an individual computational document explicitly, or changing one under `freeze: auto`, can execute its R code and requires the analytical environment and generated model artifacts described below.

This rendering path checks presentation and links. It does not independently rerun or validate the statistical analysis.

## Recompute the analyses

Use **R 4.6.1** and the versions recorded in `renv.lock`. From the repository root:

```sh
Rscript --vanilla restore.R
Rscript --vanilla run.R all
Rscript --vanilla R/run_noncyclic_sz.R
Rscript --vanilla R/run_phase_response.R
Rscript --vanilla R/sleep_diary_alignment.R
Rscript --vanilla R/manuscript_review_checks.R
Rscript --vanilla R/build_manuscript_figures.R
R_PROFILE_USER=.Rprofile quarto render --to html
```

`restore.R` installs packages into `renv/library/standalone`. The code prefers this library, so a second checkout of the parent MeLiDos project is not needed. `MELIDOS_R_LIBRARY` can specify an alternative library. On macOS, an optional OpenMP library is loaded when available; `MELIDOS_OPENMP_LIBRARY` can override its location.

The full model, validation and attribution calculations are computationally intensive and generate many GB of fitted objects. The stage driver supports `prepare`, `primary`, `chest`, `outputs`, `validation`, `attribution`, `sensitivity`, `meq`, `daylight`, `checks`, and `figures`. Saved fitted objects are reused only where the corresponding analytical checks permit it. Use `run.R figures` for display-only regeneration from existing fits and inference.

Supplement S2 has its own scenario workflow. Follow [its reproduction instructions](sensitivity/diary-msfsc/README.md), completing the matched questionnaire analysis before the diary-derived analysis and then the fixed-quartile extension. These outputs remain separate from the primary results.

After recomputation, explicitly render each changed analytical document to refresh its frozen output. Commit the updated sources, `_freeze/`, figures, tables and provenance together. `freeze: auto` detects document changes, but it does not detect every upstream data or R helper change; those changes require an explicit rerender.

## Word manuscript

With the R environment available:

```sh
R_PROFILE_USER=.Rprofile quarto render index.qmd --to docx
```

The result is `report/chronotype-manuscript.docx`. The reference document and Lua filter match the main MeLiDos manuscript style. A post-render R step applies Word layout and normalizes paragraph properties. Word review files and author-attributed redlines are local editorial artifacts and are not published on GitHub Pages.

## GitHub Pages

The [Pages workflow](.github/workflows/pages.yml) runs on pushes to `main` and on manual dispatch. It installs the pinned Quarto version, renders the website from `_freeze/`, checks local links and deploys `report/` to GitHub Pages. It does not install the research package library or refit models. GitHub Pages uses **GitHub Actions** as its publishing source.

## Inputs, provenance and attribution

Released snapshots of this repository, including the prepared inputs, analysis code, Quarto sources and retained outputs, are archived on Zenodo at [doi:10.5281/zenodo.24140946](https://doi.org/10.5281/zenodo.24140946). GitHub contains the ongoing development version; Zenodo preserves the archived releases.

The prepared inputs originate from [Zauner et al.'s MeLiDos analysis repository](https://github.com/tscnlab/ZaunerEtAl_bioRxiv_2026), primarily at commit `a5ca7e6e9968e3ca904528a2638f215dfbbd0212`. See [input-data/README.md](input-data/README.md) and [the input manifest](input-data/manifest.csv) for paths and checksums. The small demographic snapshot used for the manuscript's descriptive summary is documented by `R/manuscript_review_checks.R` and `results/provenance/manuscript_review_inputs.csv`.

The parent manuscript is *The health-relevant architecture of the everyday light exposome*, [doi:10.64898/2026.06.04.730226](https://doi.org/10.64898/2026.06.04.730226). The copied upstream data retain their [CC BY 4.0 license](input-data/upstream-metadata/LICENSE.md) and [source citation metadata](input-data/upstream-metadata/CITATION.cff). PRC coordinates and their separate source record are described in [input-data/phase_response/README.md](input-data/phase_response/README.md). Numerical PRC coordinates were digitized by the cited source authors.

Please cite this repository using [CITATION.cff](CITATION.cff), and cite the original MeLiDos study and relevant external sources when reusing their data. No journal acceptance or manuscript DOI is implied by this repository.
