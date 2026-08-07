# Scenario catalog

## Reproduction contract

Two kinds of input are included:

1. **Frozen inputs** are the exact CSVs consumed by the runs summarized in the report. Their SHA-256 hashes, row counts, exposures, and realized pooled rates are in `manifests/scenario_manifest.csv`. These files are the numerical acceptance contract.
2. **Generators** document and recreate the declared mechanisms for future batteries. The original transient generator used for the controlled-null CSVs was not retained byte-for-byte, so regenerated controlled-null samples are mechanism-equivalent but are not claimed to reproduce the old random draws. The exact old draws are therefore shipped. The corrected homogeneous-certified-null generator is deterministic and does reproduce its frozen inputs.

All controlled scenarios retain the same 120,595 frozen WALA/WAC locations and, except for the equal-exposure intervention, the same row exposures. Total original exposure is 47,566,402 loans.

## Controlled-null scenarios

| ID | Spatial truth | Row/pool effect | Exposure design | Seeds used in report | Purpose |
|---|---|---|---|---|---|
| `flat_p50_clean` | constant logit, p=.50 | none | frozen | 0--4 | baseline-probability check |
| `flat_p78_clean` | constant logit, p=.78 | none | frozen | 0--4 | clean complete null and adaptivity |
| `flat_p97_clean` | constant logit, p=.97 | none | frozen | 0--4 | saturation stress |
| `flat_p78_gauss` | constant logit, p=.78 | independent N(0,.35) logit effect per row | frozen | 0--4 | homogeneous overdispersion |
| `flat_p78_heavytail` | constant logit, p=.78 | 90/10 two-normal mixture, total variance .35, wide/core variance ratio 10 | frozen | 0--4 | heavy-tail law mismatch |
| `flat_p78_corr01` | constant logit, p=.78 | Gaussian total variance .35 with rho=.01 shared within randomized balanced 40-row blocks | frozen | 0--2 | mild dependence stress |
| `flat_p78_corr05` | constant logit, p=.78 | same, rho=.05 | frozen | 0--2 | stronger dependence stress |
| `flat_p78_gauss_equaln` | constant logit, p=.78 | independent N(0,.35) | every row n=394 | 0--2 | isolate exposure heterogeneity |
| `smooth_gauss` | nine-term polynomial outside the window, exposure-weighted quadratic in the certified interior, cosine taper between | independent N(0,.35) | frozen | 0--2 | separate legitimate low-order refinement from invented texture |

For the heavy-tail law, the core variance is `.35/(.9+.1*10)=.1842105` and the wide variance is `1.842105`; the mixture is centered and has total variance .35.

The artificial-correlation scenarios are sensitivity experiments, not issuer calibrations. Blocks were randomized to avoid deliberately aligning shared shocks with location, but they are not substitutes for real issuer identifiers.

## Corrected homogeneous certified null

`hcn_corrected_gaussian` starts from the frozen Ginnie v3 design and its nine-term polynomial truth. Inside WALA [105,175] and WAC [3.5,5.5], the truth is exactly the exposure-weighted quadratic fit over the larger WALA [80,200], WAC [3,6] region. A separable cosine taper joins the two surfaces. Each row then receives an independent N(0,.35) logit effect and exact binomial sampling. Seeds are `88000 + seed`; report comparisons use seeds 0--4.

This scenario was used for the Lindsey before/after bug measurement because it has genuine low-order structure outside the certified region while deep hierarchical surplus inside the certified region is negligible.

## Metrics

- **Gate admissions:** active, free hierarchical coefficients marked as admitted by the statistical gate.
- **Window admissions:** admitted coefficient centers inside WALA [105,175], WAC [3.5,5.5].
- **Deep window admissions:** window admissions at hierarchical depth at least 8.
- **Deep window faces:** final triangulation faces of area-depth at least 8 whose centroids lie in the window. This is a geometric closure footprint, not an FDR numerator.
- **Nonquadratic RMS:** information-weighted RMS remaining after projecting the fitted logit correction in the certified window onto a quadratic basis. It measures invented local texture after low-order structure is removed.
