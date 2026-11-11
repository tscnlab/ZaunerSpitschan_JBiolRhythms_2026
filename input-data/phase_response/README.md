# Numerical phase-response source

The CSV is an unchanged copy of the numerical coordinates deposited with Walch et al. (2025), pinned to the repository commit and blob listed in `source_manifest.json`. Its name identifies it as a prior digitization of St Hilaire et al. (2012), Figure 3A. It is not the original raw dataset. No image extraction was performed for this analysis.

The figure repeats some observations across the circular boundary. `phase_response_read_prc()` retains one 24-hour onset interval, verifies the excluded periodic copies, converts onset to midpoint, and fits the original study's two-harmonic form. The exported source-fit residuals describe variation between the digitized participant responses and the fitted curve. They are not an estimate of digitization error.

Reproduction requires only this CSV and the saved primary phase model and inference matrices. `R/run_phase_response.R` writes source points, duplicate checks, fitted coefficients, weighted indices, sensitivity results, provenance and plots. The primary model is not refitted.

The current integration uses the two extreme quartiles and centers exposure at zero light, under both published raw and drift-corrected PRCs. It integrates log10(melanopic EDI + 0.1) + 1 from the primary fitted curves. Intervals condition on the PRC and alignment. Sensitivities vary the PRC baseline and its alignment to MSFsc.
