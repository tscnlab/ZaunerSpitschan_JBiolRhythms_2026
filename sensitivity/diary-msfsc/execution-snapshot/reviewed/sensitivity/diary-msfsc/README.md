# Diary-derived corrected midsleep sensitivity

This directory holds isolated matched-questionnaire and diary-derived (dMSFsc) analyses. Neither writes to the primary project's `results/` directory. Original model producers are copied unchanged into `source-snapshot/` and each scenario's `R/`; input-data is a read-only-by-convention symbolic link to the fixed source snapshots.

From the original project root:

```sh
Rscript --vanilla sensitivity/diary-msfsc/prepare.R
```

In either `sensitivity/diary-msfsc/matched-msfsc/` or `sensitivity/diary-msfsc/dmsfsc/`, use the existing R library:

```sh
export MELIDOS_MAIN=/Users/zauner/Documents/Arbeit/12-TUM/MeLiDos/WP2.2.5_Data_Analysis/ZaunerEtAl_bioRxiv_2026
export MELIDOS_R_LIBRARY=/Users/zauner/Documents/Arbeit/12-TUM/MeLiDos/WP2.2.5_Data_Analysis/chronotype_analysis/renv/library/standalone
for stage in primary validation daylight daylight-validation prc meq sensitivity noncyclic outputs attribution
do
  Rscript --vanilla -e 'source("../scenario.R")' "$stage" > "results/logs/$stage.log" 2>&1 || break
done
Rscript --vanilla ../verify-stages.R
```

The matched primary run must finish before the dMSFsc daylight stage because the latter additionally evaluates the matched-questionnaire windows. Existing model caches are reused only when their analytical signatures match.

Then from the original project root:

```sh
Rscript --vanilla sensitivity/diary-msfsc/compare.R
R_PROFILE_USER=.Rprofile quarto render diary-msfsc-sensitivity.qmd
```

`adapters.R` generalizes two report assumptions rather than changing the fitted models: the number of significant windows is not forced to two, and unsupported daylength scenarios are calculated with the existing support mask instead of rejected. MEQ results are assembled without the original hardcoded sample-count caption/check. The noncyclic sensitivity replaces the original fixed-140 assertion with exact sample agreement against its matched cyclic baseline, retaining the original basis, endpoint, constraint, fitting and boundary checks.

Each scenario repeats the primary near-eye/chest clock/phase fits, registered inference, work/free contrasts, diagnostics, individual curves, photoperiod summaries, held-out validation, daylength and target sensitivities, daylength validation, PRC and Emens calculations, MEQ comparison, model-structure sensitivities, noncyclic sensitivity, matched placement summaries and all Shapley subsets. These outputs are not substituted into the manuscript.

Scientific definitions and limitations: dMSFsc uses recorded night counts and local-wall-clock sleep onset to final awakening. Both day types are required. No alarm information is available. The stable dMSFsc predictor uses the full diary; held-out-light prediction is therefore conditional on that diary summary. The published DLMO offsets are transferred from questionnaire MSFsc as an explicit sensitivity assumption. Selected-window average intervals are unadjusted and are kept separate from simultaneously significant windows and fixed-window comparisons.
