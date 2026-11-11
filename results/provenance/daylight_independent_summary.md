# Independent daylight verification

The preserved full-model verification script completed **37 checks**, all passing, on 22 October 2026 at 17:50 CEST. It reads the corrected cached model and does not refit models or generate coefficient draws.

Run from the secondary-analysis project root:

```sh
Rscript --vanilla results/provenance/daylight_independent_verify.R
```

The checks cover participant weighting, weighted scenario quantiles, all five training-only preprocessing splits, the numeric-by cyclic sum-to-zero basis, contrast construction, direct prediction reconstruction, exact interval counts including the circular boundary, and observed site support.

Three findings were corrected before the full-model checks passed:

- Numeric `by` handling in mgcv 1.9-4 discarded the original `sz` constraint. The local `szphoto` constructor retains the standard `sz` constructor and sets `attr(C, "always.apply")`. Its `smooth.info` method supplies the standard `tensor.possible` flag for discrete fitting. The corrected fitted term has 33 columns and its prediction matrices sum to zero across quartiles.
- Window counts used complete half-hour bins beyond the reported ten-minute interval boundaries. Counts now use the exact circular phase interval.
- Completing a factor-valued site table added artificial Q2/Q3 zero rows. Completion now explicitly includes Q1 and Q4 only. The actual data still have no site with at least five participants in both extreme groups.

The focused constructor and reduced-model check was run during the original audit and is preserved separately:

```sh
Rscript --vanilla results/provenance/daylight_independent_basis.R
```

That script reproduces the earlier `smoothCon` check and a seeded 4,000-row discrete `bam` fit. The preserved file was syntax-checked after creation; the reduced model was not refitted again during archival. Its original passing results and runtime are in `daylight_independent_basis_original_runtime.txt`.

Files:

- `daylight_independent_checks.csv`: the 37 passing checks.
- `daylight_independent_inputs.csv`: checksums of data, fitted model, fold inputs, helper code, and the executed verification script.
- `daylight_independent_folds.csv`: training and assessment participant counts.
- `daylight_independent_runtime.txt`: exact command, working directory, execution time, OpenMP path, R version, and package versions for the preserved verification run.
- `daylight_independent_original_runtime.txt`: the original full-model audit record.
- `daylight_independent_basis_original_runtime.txt`: the original constructor and reduced-model audit record.

Runtime: R 4.6.1; mgcv 1.9-4; gratia 0.11.2; dplyr 1.2.1; tidyr 1.3.2; ggplot2 4.0.3. OpenMP was preloaded from `/opt/homebrew/opt/libomp/lib/libomp.dylib`.

The audit checks the implementation and data separation. It does not replace residual diagnostics or validate a causal interpretation of photoperiod moderation.
