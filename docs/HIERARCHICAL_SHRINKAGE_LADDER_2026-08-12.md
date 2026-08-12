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
- **`eb`** (`DMESH_RESHRINK_LADDER=eb`): principled empirical Bayes. The
  prior variance is chosen by *selection-corrected marginal maximum
  likelihood* (`mle_tau`) rather than method of moments: each admitted
  surplus has calibrated sampling variance `svar = admit_sd^2` (already
  law-aware -- `admit_sd = sd_cal` carries the gate's pool term), so its
  marginal under `N(0, tau)` is `N(0, tau+svar)` truncated to the admission
  region `|draw| > athr`, and tau maximizes the summed selection-corrected
  marginal log-likelihood. The per-depth tau's are then chosen
  *hierarchically*: each depth's MMLE is shrunk toward a global MMLE by the
  precision weight `m/(m+kappa)` (`DMESH_EB_BORROW`, default 8), so sparse
  depths borrow strength instead of using the ad-hoc `4^{-dd}` decay. The
  per-vertex law-aware posterior (`tau/(1+r)`) is applied on top.
  `DMESH_EB_FREE=1` lifts the improvement-only clamp against the loop tau.
- **`ebcell`** (`DMESH_RESHRINK_LADDER=ebcell`): the same MMLE, three-level:
  `(depth x K_eff cell)` tau's borrow from the depth MMLE, which borrows
  from the global -- the evidence-geometry structure as a hierarchy level
  rather than a hard split.
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
| depth          | 1.1552  | 0.039  | pooled per-depth moment (floor-fixed default) |
| cells_rawkeff  | 1.1606  | 0.038  | hard threshold, binomial count          |
| cells          | 1.1595  | 0.033  | hard threshold, law-aware count         |
| keffcont       | 1.1534  | 0.038  | moment tau, threshold-free law-aware posterior |
| **eb**         | **1.1502** | 0.037 | **MMLE tau + hierarchical borrow (recommended)** |
| ebcell         | 1.1519  | 0.038  | MMLE tau, three-level (depth x cell)    |
| eb_free        | 1.1490  | 0.039  | eb with the improvement-only clamp lifted |
| twogroups      | 1.1671  | 0.034  | per-depth spike-and-slab                |
| twogroups_free | 1.1706  | 0.035  | slab cap lifted                         |

The ranking is the result, and it is monotone in principle -- each step
that adds statistical principle improves the held-out fit:

1. **Pool effects first beats raw counting.** `cells` (law-aware K_eff)
   improves on `cells_rawkeff` in both mean (1.1595 vs 1.1606) and spread
   (sd 0.033 vs 0.038). Deflating the effective count by the pool-effect
   variance is a measured net positive.
2. **Removing the threshold beats the threshold.** Both hard-threshold
   `cells` variants are worse than no reshrink; the threshold-free
   `keffcont` is the first ladder that beats the baselines. It wins or ties
   `depth` on 8 of 11 seeds; the small losses are seeds where the reshrink
   pass itself is mildly harmful (e.g. seed 107, where `off` is already the
   second-best fit and every reshrink arm regresses it).
3. **Marginal maximum likelihood + hierarchical borrowing beats method of
   moments.** `eb` (1.1502) improves on the moment-based `keffcont` (1.1534)
   by 0.0032 and on no-reshrink (1.1545) by 0.0043 -- the largest robust
   gain of the program. It wins on **10 of 11 seeds** (only seed 107, the
   reshrink-hostile one, loses, by 0.008), with the tightest spread of the
   net-positive arms and no catastrophic seed. The hierarchical borrow is
   the mechanism: at sparse depths where the MMLE cannot self-estimate
   (seed 103 depth 8: m=2, `tau_mle=0`) the global pull supplies a sane
   prior (`tau*=2790`), while well-populated depths with a strong
   noise signal keep their own low estimate (seed 103 depth 6: m=40,
   `tau_mle=64`, stays at 635). The three-level `ebcell` does not add over
   `eb` here (1.1519): the K_eff cell split, valuable for `cells`, is
   already subsumed by the per-vertex `tau/(1+r)` posterior once the tau is
   MMLE-chosen, so the extra level mostly adds estimation noise.
4. **Lifting the clamp trades mean for tail risk.** `eb_free` (1.1490) has a
   slightly better mean but a -0.043 blow-up on seed 104: releasing the
   improvement-only cap lets the MMLE shrink less and occasionally readmit a
   coefficient the loop had damped. `eb` (clamped) is the recommended
   default for that reason -- it captures nearly all the gain with none of
   the tail.
5. **The per-depth spike-and-slab is not paying.** `twogroups` over-shrinks
   and `twogroups_free` is worse; the EM separation is weak at shallow tiers
   (slab pi ~ 0.5, tau1 pinned at the loop tau), consistent with the
   writeup's note that the correct object is spike-and-slab per (depth x
   evidence-geometry) cell. `eb` reaches the same principled EB target
   through MMLE + hierarchy without the fragile per-tier mixture.

Mechanism check: `keffcont` marks 13-17% of admitted coefficients per run
as pool-dominated (`r > 1`, pool variance exceeds binomial), exactly the
vertices the law-aware posterior shrinks harder; `eb` keeps that posterior
and additionally chooses the tier variance by likelihood. The largest `eb`
wins land where pool-contamination is signal-scale (seed 105, +0.011 vs
off; seed 108, +0.009), confirming the gain is the principled variance
estimate, not uniform extra shrinkage.

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

`DMESH_RESHRINK_LADDER=eb` is the recommended reshrink ladder on the one-law
SMM configuration. It is the most principled of the variants -- the prior
variance is chosen by selection-corrected marginal maximum likelihood, the
tier variances are chosen hierarchically (sparse depths borrow strength
from a global fit), and the per-vertex posterior is pool-effect-aware -- and
it is the strongest measured: held-out deviance/pool 1.1502 vs 1.1545
no-reshrink and 1.1552 for the pooled-moment default, winning on 10 of 11
seeds with no catastrophic seed and the tightest spread of the net-positive
arms. Nulls and the default (non-reshrink) path are byte-identical.

The absolute magnitude is still modest (~0.004 deviance/pool); the value is
that the shrinkage coefficients are now genuinely empirical-Bayes -- chosen
by likelihood and by a hierarchy -- rather than a per-depth moment with an
ad-hoc decay borrow. `eb_free` is available for the aggressive setting but
carries tail risk (a -0.043 seed) and is not recommended as default.

The remaining target is the writeup's single spike-and-slab per (depth x
evidence-geometry) cell that unifies the admission gate and the shrinkage
ladder (bar 2.7355 on the 2020 conditional-dual release). `eb` reaches the
principled-EB variance-estimation half of that target on the current
configuration; `twogroups` is a first, not-yet-paying increment toward the
spike-and-slab half, and `ebcell` shows the cell hierarchy is already
subsumed by the MMLE + law-aware posterior for point prediction.

## Scale-exchangeable prior and honest per-pool variance (follow-on)

Two refinements motivated by the observation that K_eff measures *data
support*, not *signal*: the prior should be exchangeable over signal
coordinates (scale) and factor K_eff out through the likelihood, and the
per-pool variance `v_i` that carries all the measurement/pool information
should be honest. Harness and raw records in
`experiments/tau_floor_fix_aug12/` (`sweep_scale2.py`, `gen_pointvar.py`,
`results_scale2.json`, `analysis_scale.txt`).

- **`ebscale`** (`DMESH_RESHRINK_LADDER=ebscale`): one smooth per-scale law
  `tau(d) = tau0*rho^d` fit by MMLE over all admitted coefficients jointly,
  used as the hyperprior MEAN that each depth's own MMLE borrows toward; the
  prior never sees K_eff (it enters only through `v_i` and the `/(1+r)`
  posterior). `DMESH_EB_FREE` lifts the improvement-only clamp.
- **Honest `v_i`** via `DMESH_POINT_VAR`: replaces the shipped 4-bucket
  size-stratified pool variance (an 8.5x step at n0=512) with a smooth
  log-log interpolation in pool size (`gen_pointvar.py`), applied fleet-wide
  so it feeds both the gate null and the reshrink `r`. Pure size covariate,
  no residual conditioning, so no surface-error absorption.

Held-out deviance/pool, June 2026 SMM design, 11 seeds:

| arm             | mean   | sd     | notes                                  |
|-----------------|--------|--------|----------------------------------------|
| off             | 1.1545 | 0.039  | no reshrink                            |
| eb              | 1.1502 | 0.038  | per-depth MMLE (robust default)        |
| eb_pv           | 1.1488 | 0.040  | eb + honest v_i                        |
| ebscale         | 1.1508 | 0.036  | scale prior, clamped (tightest spread) |
| ebscale_free    | 1.1486 | 0.039  | scale prior, clamp lifted              |
| ebscale_free_pv | 1.1421 | 0.041  | scale prior + honest v_i, clamp lifted |

Reading:

1. **Confirms the exchangeability principle.** Clamped `ebscale` matches
   `eb` on mean with the tightest spread of any arm (sd 0.036): the scale
   prior's safe contribution is stability -- it stops the prior chasing
   per-depth sampling noise. Combined with the earlier result that `ebcell`
   (K_eff *in* the prior) did not pay, this is direct evidence that K_eff
   belongs in the likelihood, not the prior.
2. **Honest `v_i` is a real but modest, not-free win.** De-bucketing the
   pool variance helps the mean wherever applied (+0.0014 on `eb`, +0.0065
   on the aggressive scale prior -- more when the prior is thinner and `v_i`
   carries more weight), confirming the measurement model was costing us. It
   also lifts sd and can regress a fragile seed (s104: 1.130 -> 1.154), and
   this is only the crudest de-bucketing; per-pool predictive `E[u^2|data]`
   is the deeper version, deferred.
3. **The best mean (1.1421, +0.0125 over off) is aggressive-only.** It
   requires lifting the improvement-only clamp, which carries a tail (s104,
   1.156). Diagnosed: the admitted-only selection-corrected MMLE over-reads
   a pure-spike deep tier as signal (`tau_mle` 6428/13191 at depths 8-9
   where held-out wants the floor), and only the clamp catches it. This is
   the signature of fitting `tau` on selected coefficients and re-motivates
   the full-ensemble spike-and-slab (which would estimate the null fraction
   pi directly and not need the clamp).

### Recommended default and next direction

`eb` remains the robust default; `ebscale` is available where lower
cross-seed variance is worth more than the last hundredth of mean, and the
aggressive `ebscale_free_pv` is the measured ceiling (1.1421) but not yet
bankable. The clamp-vs-free tension should not be resolved by tuning the
clamp. The planned next step is an **EB tau multiplier**: replace the hard
`min(tau_reshrink, tau_loop)` with a continuous factor `c` (per tier)
applied to the loop tau, estimated by marginal likelihood and regularized
toward 1, so the data chooses how far to move off the in-loop
regularization and the s104-style over-read is damped by the pull toward 1
rather than caught by a kink. That, together with fitting the prior on the
full candidate ensemble (refusals included) rather than admitted-only,
is the route to capturing the aggressive gains safely.
