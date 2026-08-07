# Pooled-permutation gate: frozen null-battery results

Date: 2026-08-01

## Test set

The battery uses 34 frozen complete-null datasets on the real Ginnie Mae
WALA/WAC locations and exposure design:

- clean binomial nulls at p = .50, .78, and .97 (five seeds each);
- independent Gaussian pool effects (five seeds);
- heavy-tailed pool effects (five seeds);
- randomized block dependence at rho = .01 and .05 (three seeds each);
- Gaussian effects with equal exposure (three seeds).

All runs use q = .10, the local held-out split, three fit stages, parent variance
scale .25, and the current exact-fit C++ branch.

The experimental permutation implementation uses 32 deterministically seeded
permutation draws per candidate, pools the null draws across candidates by star
size, and applies direct BH. It is not yet the memo's deterministic
double-saddlepoint implementation.

## Three-round battery

| Rule | Mean false admissions | Runs with any | Maximum |
|---|---:|---:|---:|
| Shipped gate | 17.50 | 34/34 | 158 |
| Shipped gate + K_eff >= 4 | 29.24 | 34/34 | 224 |
| Permutation BH, K_eff >= .5 | 0.029 | 1/34 | 1 |
| Permutation BH, K_eff >= 4 | 0.000 | 0/34 | 0 |

The support floor alone does not explain the improvement: raising K_eff without
the permutation calibration did not control the battery and increased the
overall mean admissions in this adaptive implementation. With the original
K_eff = .5 threshold, the permutation rule produced one false admitted
coefficient in 34 runs. Raising the threshold to 4 removed that final
admission.

### Scenario means

| Scenario | Shipped | Permutation BH, K_eff >= .5 | Permutation BH, K_eff >= 4 |
|---|---:|---:|---:|
| Clean p=.50 | 5.20 | 0.00 | 0.00 |
| Clean p=.78 | 37.00 | 0.00 | 0.00 |
| Clean p=.97 | 6.80 | 0.00 | 0.00 |
| Gaussian pool | 16.00 | 0.20 | 0.00 |
| Heavy-tail pool | 18.40 | 0.00 | 0.00 |
| Correlated ρ=.01 | 20.67 | 0.00 | 0.00 |
| Correlated ρ=.05 | 20.33 | 0.00 | 0.00 |
| Gaussian equal n | 18.33 | 0.00 | 0.00 |

## Full adaptation, up to 80 rounds per stage

| Rule | Mean false admissions | Runs with any | Maximum |
|---|---:|---:|---:|
| Shipped gate | 121.94 | 34/34 | 1001 |
| Permutation BH, K_eff >= .5 | 0.029 | 1/34 | 1 |
| Permutation BH, K_eff >= 4 | 0.000 | 0/34 | 0 |
| Permutation theoretical-null lfdr, K_eff >= 4 | 0.000 | 0/34 | 0 |

The single K_eff=.5 false discovery remained a single final coefficient after
the full adaptive run; this battery did not exhibit descendant amplification.
That is encouraging but is not a formal branchwise error-spending result.

## Assessment

The pooled-permutation calibration passes this frozen complete-null battery
decisively. The main open problem is now power rather than null size: on the
real Ginnie held-out fit, the experimental gate previously retained too little
structure and scored worse than the shipped estimator.

The direct-BH form is the safer experimental baseline. The theoretical-null
lfdr form also made no false admissions here, but its real-data behavior and
the earlier isolated heavy-tail instability make it less suitable as the
default until investigated further.

## Important limitations

- 34 null datasets are enough to reveal the large shipped-gate failure but only
  give coarse resolution for a target familywise error probability near .10.
- The current implementation uses pooled deterministic pseudo-random draws,
  not the proposed double-saddlepoint tail calculation.
- The battery tests complete-null behavior; it does not measure power or
  selection quality when genuine spatial structure is present.
- Although the 80-round experiment showed no false-branch cascade, branchwise
  error spending remains unresolved theoretically.
