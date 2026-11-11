Current PRC and daylength verification:
Rscript --vanilla results/provenance/revision4_verify_analysis.R
This independently reconstructs observed daily medians, equal participant weights, all 120 indices, all 60 contrasts and their full coefficient covariance; checks daylength targets, support and linear scaling; and preserves protected model caches.
The phase_response_independent_source.R verification remains applicable to the unchanged published coordinates and fitted PRC.
Earlier integration/result/claim audits are archived in revision4-original/results/provenance and describe the superseded four-group, uncentered/clipped PRC analysis.
The source-fit audit retains an adjusted-R-squared discrepancy: conventional 18-row OLS gives 0.454 versus the published 0.55. Published extrema and amplitude reproduce; we do not claim to reproduce its adjusted R-squared.
