# Session notes, 2026-08-05/06 (overnight)

Writeup section: sec:aug6 in the main document (sec_aug6_session.tex). Full
narrative and every number lives there; this folder carries the artifacts.

- law2020.json / law2020_str.txt: converged 2020-regime SMM pool law
  (per-stratum pi, tail var, total). Standing rule: law fixed point per
  regime and per model (post-boosting refit required).
- relax_census_final.csv: one-coordinate relax duals (law-marginal currency,
  tau=1), all 270 scoreable vertices of the Jun->Jul fit; sharp/diffuse
  decomposition by held nats per ring pool.
- relax_census_conditional.csv: same duals in the conditional (binomial)
  currency; input to the hierarchical-tau conditional slab (66v release,
  held dev 2.7355 vs 2.7936 shipped).
- release_h.json: the 5-vertex local release heights (knob -1.09 -> -4.42).
- diag_jj.csv: region event log (DMESH_BRANCH_DIAG with the new
  DMESH_DIAG_REGION env patch in core/pipeline.cpp).
- hmt.py: HMT/two-groups trust prototype (audit overlay; post-hoc surface
  editing refuted, see writeup).
- pngs: key session figures (more in the writeup directory).

Binary patches in core/pipeline.cpp this session: DMESH_DIAG_REGION
(env-configurable diagnostic region box) and DMESH_MIN_KEFF (candidate
support floor; refuted, default off, retained for reproducibility).

Bug lead: TAU_FLOOR bypass at depths 8-9 (lambda_v = 1e8 wholesale vs the
floor-implied cap of 4; d10 shows the floor working). Source read queued.

Source read (2026-08-07): RESOLVED — defect, two-pass interaction, fixed.
The DMESH_FINAL_RESHRINK block sets tau_floor_sq = 1e-8 for its own no-floor
sweep and never restores it. The subsequent retirement pass
(retire_zero_support_at_fixed_topology), whenever it retires anything, calls
recompute_tau_sq: a fresh in-loop DL solve over ALREADY-SHRUNK delta_pooled
(chi ~ 0 at strongly shrunk depths) under the gutted floor -> tau = 0 ->
max(0, 1e-8) -> lambda = 1/1e-8 = 1e8 exactly, wholesale per depth, in the
post-retirement sweep and final refit. It also clobbers the reshrink taus.
Reproduced locally on seed7-lindsey + DMESH_FINAL_RESHRINK=1 (v5145, d10,
keff 2.3, lambda 1e8). Fix validated here (core/pipeline.cpp): save/restore
tau_floor_sq around the reshrink sweep. Result: weld gone, lambda = 4 (the
floor cap, as documented), held-out unchanged (6.455 vs 6.454 dev/pool).
Adopted implementation (rebase note, 2026-08-13): the Aug 12 session fixed
the same defect at the root (a8c18d2) — the reshrink pass now respects the
configured floor outright, with DMESH_RESHRINK_NO_FLOOR=1 preserving the
historical no-floor behavior — superseding the save/restore variant, which
was dropped in the rebase. Same mechanism, same result.

Related edge case (empty-star vertices, same session): support-starved
(contained) vertices — stars over near-empty space, keff < 0.5 or
fit_points < 10 — are frozen at their prior by update_height but still vote
in both tau moments (recompute_tau_sq and the reshrink collection) with
stale admission sigma and noise surpluses: they deflate tau where their
deltas are contained (seed7 d10: 2 of 3 voters are keff~0.01, chi=0.3 vs
m-1=2 -> tau=0 -> the one real workhorse welded at lambda=4) and inflate it
where they are wild (d4/d7: tau 2.5x/1.5x high). New env patch
DMESH_TAU_EXCLUDE_CONTAINED=1 (default OFF, reproducibility) removes their
vote; contained vertices stay fully pooled to the prior. On seed7-lindsey
the exclusion frees the d10 workhorse (real surplus, sane lambda 3e-4) but
reads slightly WORSE held-out (6.478 vs 6.455; topology feedback, 1873 vs
2196 vertices) — fresh-seed confirmation required before adopting.
