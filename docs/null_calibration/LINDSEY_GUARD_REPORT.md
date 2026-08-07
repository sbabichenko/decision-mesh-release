# Exact-score gate continuation: invalid Lindsey fit guard

## Finding

The catastrophic B3 branch on corrected homogeneous-null seed 4 was caused in large part by a numerical failure in the Lindsey Poisson histogram fit, not by the exact coordinate solver.

At stage 1, round 14:

- candidate count: 988;
- score mean / sd: 0.0588 / 1.3157;
- the Poisson IRLS did not converge within 100 iterations;
- the fitted histogram density integrated to about 4322.6 rather than 1;
- central matching had nonnegative quadratic curvature;
- the code reset only `(delta0, sigma0, pi0)` to `(0,1,1)` but still used the invalid fitted density;
- resulting local FDRs were approximately 0.0000--0.0018, so all 988 candidates were admitted.

## Code change

`lindsey_lfdr` now fails closed unless:

1. the Poisson IRLS converges;
2. all coefficients and fitted bin densities are finite;
3. the fitted density normalizes to within 2% of one;
4. the central log-density fit has strictly negative quadratic curvature.

When any check fails, `pipeline.cpp` uses direct two-sided normal-score BH for that candidate family. The invalid local-density fit is never allowed to authorize admissions. Logs identify this as `BH-fallback`.

## Regression results

### Corrected homogeneous certified null, seeds 0--4

| Metric | Original B3 | Guarded B3 |
|---|---:|---:|
| Mean global gate admissions | 644.2 | 233.6 |
| Mean window admissions | 53.4 | 14.2 |
| Mean depth>=8 window admissions | 7.0 | 0.0 |
| Mean deep window faces | 268.8 | 101.0 |
| Mean nonquadratic fitted RMS | 13.62 | 6.34 centi-logits |

The known runaway seed 4 changed from 1,299 global admissions, 141 window admissions, 34 deep window admissions, 754 deep faces, and 28.42 centi-logits RMS to 174, 2, 0, 41, and 1.48, respectively.

### Other checks

- Clean simulated complete null: 81 active vertices, zero admissions, zero FDP; passed.
- Frozen v3: no invalid Lindsey fallback occurred; 80 deep window faces, so the guard does not solve the remaining calibration problem.
- Frozen v4: 2 fallback rounds, 3 depth>=8 window admissions but 287 deep window faces; still far outside the 3--7 face target.

## Decision

Accept the guard as a correctness fix. Do **not** accept guarded B3 as certified. It removes a fail-open numerical pathology and most of the seed-3/4 runaway behavior, but adaptive repeated gating still creates too much geometric footprint and fitted texture.

The next architectural step remains a fixed candidate hierarchy with independent grouped assessment and one final structured gate. Re-fitting an empirical null inside every adaptive round is too fragile to carry the error-control claim.


## Packaged reproduction scenario

The exact five corrected homogeneous-certified-null CSVs used in this comparison are now in `experiments/corrected_homogeneous_certified_null/frozen_inputs/`. The portable generator, source instrument, truth coefficients, before-guard source, guarded source, and both summary tables are included. This makes the 644.2-to-233.6 comparison a fully named before/after regression rather than an undocumented external run.
