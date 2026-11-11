# Prepared MeLiDos inputs

This is a fixed local copy of the inputs used by the chronotype secondary analysis. The source is [ZaunerEtAl_bioRxiv_2026](https://github.com/tscnlab/ZaunerEtAl_bioRxiv_2026), commit [`a5ca7e6e9968e3ca904528a2638f215dfbbd0212`](https://github.com/tscnlab/ZaunerEtAl_bioRxiv_2026/tree/a5ca7e6e9968e3ca904528a2638f215dfbbd0212), copied on 22 October 2026. The files are byte-for-byte copies, not newly calculated versions. The source checkout was clean when copied.

The 14 source files total 4,933,712 bytes. `manifest.csv` records the source repository and commit, original relative path, local relative path, size, MD5 and SHA-256 for each file. `R/data.R` checks the copied-file MD5 values before import. Data preparation reads these local files and has no dependency on the source checkout.

## Contents

| File or folder | Role |
|---|---|
| `model_data/normalized_inputs/chronotype.rds` | Participant-level MSFsc, MEQ scores and recorded MEQ categories. |
| `model_data/normalized_inputs/sleepdiaries.rds` | Sleep diaries, waking dates and reported work/free-day assessments. |
| `model_data/H02/main__glasses__all_available.rds` | Prepared near-eye half-hour observations, arithmetic-mean melanopic EDI, eligibility, timestamps and temporal provenance. |
| `model_data/H02/main__chest__all_available.rds` | Corresponding chest-level observations. |
| `model_data/base/metrics_*_30_minute_context.rds` | Matching sensor-bin context, including photoperiod, civil dawn/dusk, solar noon and site-level geographic/solar provenance. |
| `model_data/temporal_provenance/wall_outcome_links.rds` | Original elapsed-time links for auditing wall-clock bins and daylight-saving transitions. |
| `upstream-source/` | The source project's daily-pattern data preparation, model and dominance scripts, retained for methodological traceability. These scripts are not executed by this analysis. |
| `upstream-metadata/` | Unmodified upstream license, citation metadata, README and complete package lockfile. |

The source repository calls its daily-pattern analysis **H02**. That identifier remains in imported filenames so each file can be traced to its origin. This secondary project supplies its own analysis code and report.

## Attribution and licensing

The upstream `LICENSE.md` and `CITATION.cff` declare Creative Commons Attribution 4.0 International (CC BY 4.0). Their exact copies are retained under `upstream-metadata/`. Cite the MeLiDos data sources and the Zauner et al. repository/preprint when reusing these prepared data. The copied citation file preserves the complete source authorship and attribution metadata. The snapshot accompanies this reproducible analysis.

## Software snapshot

`upstream-metadata/renv.lock` is the complete original lockfile. The project-root `renv.lock` retains the existing records for 120 analysis packages and transitive dependencies. `package-manifest.csv` lists those records and separates direct runtime requirements from dependencies. Package versions were not upgraded. `restore.R` restores this environment explicitly into the standalone project library.

The data snapshot is versioned separately from the manuscript. The source commit identifies the exact prepared files used here even when a newer preprint or repository commit becomes available.
