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
