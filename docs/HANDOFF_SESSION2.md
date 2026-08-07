# Session-2 handoff: phase B (gate on exact gains + null recalibration)

## Where things stand (verified, pins in EXPERIMENT_NOTES.md)
- Your session-1 clamp discovery was verified in the legacy tree and
  the clampless line was adopted as mainline. Your lifecycle fixes
  were necessary but not sufficient: a 20-split holdout battery found
  split-dependent runaway eruptions (15% rate) that no amplitude
  clamp, cap, or scoped exact commit eliminated.
- Root cause chain (full forensics in notes): unbounded IRLS working
  responses at saturated probabilities (steps of hundreds of logits
  injected into mesh.values) seeded eruptions; gate admissions then
  locked the noise into topology (eruption vertices 75-90% admitted).
- PHASE A (in this tree, DMESH_PL_SOLVER=1): the stage height solver
  is now coordinate ascent on the exact tempered likelihood + the EB
  prior (mu_prior, lambda = 1/tau(depth)) on hierarchical columns
  with live per-pool eta. Bounded-above objective, monotone ascent:
  eruption impossible by construction. Battery: 0/20, mean 2.9444,
  sd 0.0067, worst 2.9578. The mesh now beats the same-protocol GAM
  on EVERY split (GAM 2.9808 +- 0.0045). MAP for one-sided stars is
  native (the prior term in the Newton). Gate untouched: phase-A
  contract.

## Your mission: phase B
The gate still prices admission/refinement in quadratic-surrogate
gains from the refresh snapshot: two currencies, one economy. This
mismatch caused the spike field, the admitted checkerboards, and the
edge undervaluation. Phase B:
1. GAIN CURRENCY: admission + refinement + heap scores become exact
   tempered log-likelihood improvements. pl_solver's Newton computes
   the objective change as a byproduct; the exact_score es.W
   machinery scores proposals. Mostly rerouting, not new math.
2. NULL RECALIBRATION (your session-5 plan, now scheduled): the lfdr
   null gain distribution is a property of solver + currency
   together. Rerun the certified-null battery (v3, v4, cos g0,
   chirp, confound) under the exact solver with exact gains, refit
   the null distribution, re-derive thresholds.
3. EB TAU: the truncation correction in the tau ladder follows the
   admission threshold; move them together.
Do 1-3 as ONE event so calibration is derived once, against the
final system.

## Acceptance criteria (all four, no exceptions)
- v3 deep-in-window ~ 2; v4 at pin (3-7)
- t_split: 0/20 eruptions on the 20-seed battery (script below),
  sane sd within ~2x of 0.0067
- base WALA 0-2 (6% band) fitted <= -6.9 vs law-marginal dot -7.95
- NLL all within 0.01 of 2.9347 on the original split
The stale-vertex census (stale_census_cl_def.csv) and the 404/202
eruption seeds are regression sets: re-run both.

## Tools in this tree
- t_split.py: the 20-seed battery (mesh run + mixture NLL), ~4 min.
- eval_harness.py: NLL trio + law-marginal base profile (the only
  approved instruments; empirical logits are banned program-wide).
- Env reference: DMESH_PL_SOLVER=1 (phase A solver),
  DMESH_HIER_FIT=2, DMESH_HEIGHT_LIMIT_LOGIT=12 (safety net, likely
  removable: test), DMESH_IRLS_STEP_CAP=4 (same),   DMESH_HIER_MAP (subsumed by pl_solver: retained for A/B),
  DMESH_CONTAIN_HEIGHT=1 / DMESH_PL_SOLVER=0 etc restore legacy for
  bit-exact guards (see notes for the guard table).
- SCAFFOLDING-REMOVAL TEST (queued, do early): with pl_solver on,
  the caps/clip/sat-exact should be inert. Verify battery holds with
  each off, then delete dead mechanisms rather than shipping them.

## Files
EXPERIMENT_NOTES.md carries the complete session history including
every falsified hypothesis and pin. Battery results:
bootstrap_battery*.csv. Key figures in /mnt/user-data/outputs:
eruption_schematic.png, surrogate_vs_exact.png, instrument_profiles,
regional_map_default, regime_straddle, young_region_sixpanel.
