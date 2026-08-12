# Self-child formulation A/B - 2026-08-12

Mode: `DMESH_SELF_CHILD=1` (sequential frozen-increment architecture; see
`docs/SELF_CHILD_FORMULATION_2026-08-12.md`). All runs use the shipped
`run_ginnie_benchmark.sh` lindsey configuration (PL solver, B3 score, quad
prior, local split, q=0.10) with only the flag toggled. Verdicts follow the
repo's status convention: everything below is MEASURED.

## Real Ginnie cross-section (April proxy, 120,595 pools)

Ten matched local splits, seeds 101-110, held-out mean deviance/pool
(`ginnie_10seed_heldout_dev.csv`):

- baseline mean 6.638, self-child mean 6.716;
- paired diff +0.0775 (sd 0.054, se 0.017): self-child WORSE in 9/10 splits.

For scale: the pooled-permutation gate lost the matched-split audit by ~0.31
in 10/10 and was not adopted. Self-child v1 loses by ~0.08 in 9/10. The
joint solve's continuous re-fitting is worth about 1.2% held-out deviance
on this cross-section over frozen sequential increments, in this v1 config.

Seed 7 single run: 6.455 vs 6.517; self-child wall time 7.6s vs 30.1s
(no per-round full-mesh sweeps, no 30-cycle final refit). Final meshes are
similar size (571 vs 604 active); self-child admits gated re-adjustments
down to depth 3 (base-level repricing through the gate, the intended
mechanism) but reaches less deep (no depth-10 admissions vs 3 baseline).

## Simulator (sigma_u = 0.30, amplitude default, 10 seeds, local split)

`sim_10seed_sigma030.csv`. Baseline exhibits its known topology-collapse
mode on 2/10 seeds (held-out deviance 30.1 and 47.4 with starved admission
counts 104 and 35); self-child collapses on NONE (range 15.2-20.6). On the
8 non-collapsed seeds the two are statistically tied (18.06 vs 18.18).
Empirical FDP stays far below q in both modes.

## Complete null (flat truth, no pool effects, 5 seeds)

Zero admissions in both modes, all seeds: the gate stays closed under the
enlarged (self-delta) candidate family. The lindsey fail-closed regression
test passes.

## Reading

The formulation delivers its structural claims at a measured price:

- bought: no topology-collapse mode (sequential bounded commits cannot
  erupt); ~4x wall-time; stable decision identities (a coefficient is tested
  once against a fixed offset and never silently re-fit, so admission
  certificates never mutate and coarse repricing is FDR-charged and logged);
  the coupled-penalty machinery dissolves (prior anchors are constants).
- paid: +0.078 +/- 0.017 held-out deviance/pool on the real cross-section.

## v1 caveats that plausibly contain the gap

(1) tau tiers are still keyed by geometric depth, not by round, so delta
priors inherit cumulative-surplus scales; (2) the candidate null retains the
0.25*sp parent-variance term although parents are frozen (over-wide null,
under-admission at depth); (3) self-child stalls out of rounds earlier
(stall sweeps are no-ops by design; only the IRLS relinearization can
unlock); (4) late-discovered coarse structure must be expressed as many
small deltas (the diffuse-class channel). None of these are structural to
the formulation; all four are measurable follow-ups.
