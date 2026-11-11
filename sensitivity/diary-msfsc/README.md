# Diary-derived corrected midsleep sensitivity

This directory holds isolated matched-questionnaire and diary-derived (dMSFsc) analyses. Neither writes to the primary project's `results/` directory. Original model producers are copied unchanged into `source-snapshot/` and each scenario's `R/`; input-data is a read-only-by-convention symbolic link to the fixed source snapshots.

From the original project root:

```sh
Rscript --vanilla sensitivity/diary-msfsc/prepare.R
```

From the project root, select the restored R library, then enter either `sensitivity/diary-msfsc/matched-msfsc/` or `sensitivity/diary-msfsc/dmsfsc/`:

```sh
export MELIDOS_R_LIBRARY="$(pwd)/renv/library/standalone"
cd sensitivity/diary-msfsc/matched-msfsc # repeat later in ../dmsfsc
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

## Alignment with questionnaire groups retained

The focused `dmsfsc-fixed-groups/` extension retains the matched questionnaire quartiles and changes phase alignment to dMSFsc. Its light records, site centers, day-type weights and model formulas match the existing analyses. The shared quartile assignment step is overridden locally to use the stored questionnaire value for classification while retaining diary timing. The frozen model producers are unchanged.

This extension repeats the near-eye primary phase and expanded daylength models, their inference and diagnostics, and all three photoperiod target specifications. It does not repeat clock-time, validation, attribution or PRC analyses. The comparison uses the original fixed phase windows and preserves 95% pointwise daylength bands.

From the project root, after the two complete matched scenarios exist:

```sh
Rscript --vanilla sensitivity/diary-msfsc/prepare-fixed-groups.R
```

Then in `sensitivity/diary-msfsc/dmsfsc-fixed-groups/`, with the same library environment variable as above:

```sh
Rscript --vanilla -e 'source("../fixed-groups.R")' primary > results/logs/primary.log 2>&1
Rscript --vanilla -e 'source("../fixed-groups.R")' daylight > results/logs/daylight.log 2>&1
```

Finally from the project root:

```sh
Rscript --vanilla sensitivity/diary-msfsc/compare-fixed-groups.R
R_PROFILE_USER=.Rprofile quarto render diary-msfsc-sensitivity.qmd
```

Outputs named `fixed_groups_*` compare questionnaire timing and groups, diary timing with questionnaire groups, and diary timing and groups on the same 135 participants. The matched questionnaire cutoffs are refitted on those 135 people; three retained participants therefore differ from their original 140-person questionnaire quartile. The intermediate scenario preserves the matched reference membership exactly. Comparing the first two scenarios isolates alignment; comparing the latter two isolates subsequent reclassification. This is a sequential specification comparison, not a formal test of the difference between model effects.
