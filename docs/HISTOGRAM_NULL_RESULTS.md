# Deterministic pool-law null experiment

## Implementation

- Candidate-specific score null from the exact candidate column and frozen likelihood temper.
- Four pool-size strata.
- Each stratum uses a two-component zero-mean Gaussian scale-mixture law with configurable total variance, wide-component probability, and variance ratio.
- The most influential mixture rows are enumerated exactly (default top 8); the remaining contribution is moment-matched Gaussian.
- Two-sided candidate p-values are converted to signed Gaussian-equivalent scores for diagnostics.
- Admissions use candidatewise Benjamini-Hochberg directly; the old mean-lfdr prefix is not applied a second time.
- No Monte Carlo is used.

## Certified-window counts

Acceptance targets from the null-instrument README: v3 about 2; v4 3--7.

| Configuration | v3 deep false faces | v4 deep false faces |
|---|---:|---:|
| B3 sandwich + gate variance 0.35 | 80 | unfinished in prior run |
| Histogram null, Gaussian pool law | 80 | not run |
| Histogram null, measured two-component shape | 40 | 78 |
| Histogram null, stronger/full tail law | 35 | not run |

The measured histogram law cuts v3 false refinement roughly in half, but does not certify the exact gate. On v4 it remains far above the target.

## Real-data pin

The full-tail histogram gate produced:

- held-out marginal NLL: 2.9620;
- WALA 0--2 fitted logit: -5.68 (law-marginal reference -7.95).

This is materially worse than the best B3 predictive/boundary line (about 2.928 NLL and -6.99 at the base).

## Verdict

Pool-law tails explain a meaningful part of the false-refinement problem, but not most of the remaining v4 failure. The deterministic candidatewise histogram null is therefore a useful diagnostic and possible component of a later gate, not a complete replacement yet. Remaining likely causes are issuer/shared dependence, selection-conditioned candidate families, and local fit-state bias not represented by independent pool-effect histograms.
