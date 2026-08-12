# Hierarchical shrinkage: evidence-geometry ladder variants (2026-08-12)

Follows the tau-floor fix (`docs/TAU_FLOOR_RESHRINK_FIX_2026-08-12.md`) and
attacks the structural defects the August 6 EB conceptual audit named
(writeup sec:aug6-ladder): depth is the wrong exchangeability unit, and the
learning signal is corrupted by pool effects the working-binomial noise
denominator ignores. All variants live in the final post-selection
re-shrinkage pass (`DMESH_FINAL_RESHRINK=1`) and are selected by
`DMESH_RESHRINK_LADDER`; the default (unset) is the pooled per-depth
estimator with the restored floor, unchanged.

## Variants

- **`cells`** (`DMESH_RESHRINK_LADDER=cells`): split each depth at a K_eff
  threshold (`DMESH_TAU_CELL_KEFF`, default 1.0) and run the
  truncation-corrected moment separately in each cell, so a noise-dominated
  starved cell cannot set the workhorse cell's tau (EB audit defect 1).
- **Law-aware K_eff** (default inside `cells`): the split thresholds on
  `keff_law = keff / (1 + r)`, `r = Var_pool / Var_binom`, the pool-effect
  variance inflation the admission gate already prices (one-law `xT2xs`
  under `DMESH_ONELAW`, flat `xT2x` model otherwise). A tier of a few
  mega-pools -- each a single draw of its pool effect -- no longer reads as
  high surface support (EB audit defect 2). `DMESH_TAU_CELL_RAWKEFF=1`
  restores the binomial count to isolate the deflation.
- **`keffcont`** (`DMESH_RESHRINK_LADDER=keffcont`): the threshold-free
  continuous limit of `cells`. No K_eff cut; each vertex takes the honest
  per-depth tier tau divided by its own contamination, `tau_i = tau_depth /
  (1 + r_i)`. Pool-dominated vertices shrink smoothly toward the floor,
  clean surface witnesses keep the full tier variance.
- **`twogroups`** (`DMESH_RESHRINK_LADDER=twogroups`): per-depth
  spike-and-slab EM on the raw surpluses (spike at the floor, slab from the
  responsibility-weighted truncation-corrected moment, selection through the
  admission normalizer); each coefficient gets the posterior-mixture prior
  variance `gamma*tau1 + (1-gamma)*tau0` (EB audit defect 3, first
  increment). `DMESH_RESHRINK_2G_FREE=1` lifts the improvement-only cap on
  the slab.

Per-vertex prior variances are carried by a new `Vertex::tau_override_sq`
consumed by `compute_eb_params` in the fixed-topology re-solve. All variants
respect the restored tau floor.

## A/B measurement

Harness and raw results: `experiments/tau_floor_fix_aug12/` (`sweep2.py`,
`sweep_kc.py`, `analyze2.py`, `results_ladder.json`, `analysis_ladder.txt`).
`off` = no reshrink; `depth` = the current pooled estimator (post-floor
fix). Held-out mean deviance/pool, June 2026 SMM design, 11 seeds
(7, 101-110):

| arm            | mean    | sd     | notes                                   |
|----------------|---------|--------|-----------------------------------------|
| off            | 1.1545  | 0.039  | no reshrink                             |
| depth          | 1.1552  | 0.039  | pooled per-depth (floor-fixed default)  |
| cells_rawkeff  | 1.1606  | 0.038  | hard threshold, binomial count          |
| cells          | 1.1595  | 0.033  | hard threshold, law-aware count         |
| **keffcont**   | **1.1534** | 0.038 | **threshold-free, per-vertex**        |
| twogroups      | 1.1671  | 0.034  | per-depth spike-and-slab                |
| twogroups_free | 1.1706  | 0.035  | slab cap lifted                         |

The ranking is the result, and it is monotone in principle:

1. **Pool effects first beats raw counting.** `cells` (law-aware K_eff)
   improves on `cells_rawkeff` in both mean (1.1595 vs 1.1606) and spread
   (sd 0.033 vs 0.038). Deflating the effective count by the pool-effect
   variance is a measured net positive.
2. **Removing the threshold beats the threshold.** Both hard-threshold
   `cells` variants are worse than no reshrink; the threshold-free
   `keffcont` is the only ladder that beats the baselines, improving on
   `off` by 0.0011 and on `depth` by 0.0018. It wins or ties `depth` on
   8 of 11 seeds; the three small losses (<=0.003) are all seeds where the
   reshrink pass itself is mildly harmful (e.g. seed 107, where `off` is
   already the second-best fit and every reshrink arm regresses it).
3. **The per-depth spike-and-slab is not yet paying.** `twogroups`
   over-shrinks relative to the continuous law-aware prior, and lifting the
   slab cap (`twogroups_free`) makes it worse. The EM separation is weak at
   shallow tiers (slab pi ~ 0.5, tau1 pinned at the loop tau), consistent
   with the writeup's note that the correct object is a single
   spike-and-slab per (depth x evidence-geometry) cell rather than per
   depth. Left in as an experimental flag, not recommended.

Mechanism check: `keffcont` marks 13-17% of admitted coefficients per run
as pool-dominated (`r > 1`, pool variance exceeds binomial), and those are
exactly the vertices it shrinks harder. The largest wins land where that
set is real signal-scale contamination (seed 105, +0.011 vs off; seed 108,
+0.007), confirming the gain comes from the pool-effect deflation rather
than uniform extra shrinkage.

### Null and April-proxy behavior

Shipped flat true-null (`aprilnull`) is unchanged across all ladder
variants except `twogroups_free` (+0.002 on one split); `keffcont`,
`cells`, and `twogroups` equal the baseline. On the April proxy (flat
`GATE_EXTRA_VAR` model, constructed attrition response) the deflation uses
the homogeneous `xT2x` term so per-vertex variation is small; `keffcont` is
neutral-to-slightly-worse (+0.001 mean) there. The gains are specific to
the one-law SMM regime with a real per-pool variance object, which is the
deployment target.

## Recommendation

`DMESH_RESHRINK_LADDER=keffcont` is the recommended reshrink ladder on the
one-law SMM configuration: it is the most principled (per-vertex,
threshold-free, pool-effect-aware), the only variant that beats no-reshrink,
and it leaves the null and default paths intact. The magnitude is modest
(~0.002 deviance/pool over the pooled estimator); the structural value is
retiring the depth-only exchangeability unit and the arbitrary K_eff cut
in one continuous law-aware object. The per-cell spike-and-slab remains the
writeup's target design (bar 2.7355 on the 2020 conditional-dual release);
`twogroups` here is a first per-depth increment toward it, not the finished
object.
