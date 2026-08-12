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

## Addendum (same day): gap diagnosis and the refresh-cadence fix

The +0.0775 v1 deficit decomposes into three named parts (all MEASURED,
seed 7 forensics + 10-seed confirmation):

1. **Stale scoring linearization (~half the gap).** The baseline's per-round
   joint sweeps double as relinearization for the candidate scorer; freezing
   heights silently froze the scoring epoch. `DMESH_REFRESH_EVERY=1`
   (relinearize every round, heights untouched) fixes it. Stall patience
   alone (`DMESH_STALL_LIMIT`) recovers nothing: scores cannot change
   without relinearization. 10-seed paired result with rf1:
   all-pools +0.0369 (sd 0.071, se 0.023), vs +0.0775 (se 0.017) at v1
   (`ginnie_10seed_rf1_heldout_dev.csv`). Null battery at rf1: still zero
   admissions.

2. **WALA-0 artifact chasing.** The seed-7 gap concentrated in the issuance
   boundary at high coupon; 1,651 WALA-0 pools alone carry +0.81/pool - the
   class the writeup's own audit excludes as proxy-unreliable. Under the
   writeup's exclusion protocol the 10-seed rf1 gap is +0.0305 (se 0.011,
   8/10 worse) - smaller and much lower-variance than the all-pools metric.

3. **Renegotiation debt (the residual, concentrated in steep-gradient deep
   refinement).** In Gauss-Seidel, child evidence flows upward: parents
   track E[.|current model] continuously, so coefficient trajectories are
   martingale-like and deep candidates are priced against a locally
   re-optimized backdrop. Frozen increments accumulate the upward-flow
   correction as predictable drift in sub-threshold pieces the gate never
   books; the corner (a ~1 logit / 4 WALA-month ramp under huge per-pool
   information) is where the debt compounds: SC final corner mesh has
   d9:7-18/d10:0 vs baseline d9:34/d10:3. Refuted alternatives, each by
   direct experiment: EB over-shrinkage (shrink factors 0.93-0.99),
   tau-floor sensitivity (zero effect), stall budget alone (zero effect),
   parent-variance null width (helps only under stale refresh; hurts with
   rf1), quadratic prior target (removing it hurts).

   **Naive repair refuted:** `DMESH_SC_RENEG` thaws admitted vertices'
   parents (1) or free one-ring (2) into the round polish, un-gated, prior
   centered at no-change. Both grades WORSEN seed 7 (6.595 / 6.533 vs 6.486
   plain rf1): a parent renegotiating without its full co-support family in
   the solve fixes the admitted child's room and silently damages the other
   room - the constrained-closure compromise reproduced from the opposite
   direction. The flag is retained as a measured-negative experiment.

   The open repair is therefore closure-complete or gated: group-move
   candidates (score the joint re-solve of a parent's full free family as
   ONE gated decision - the aug6 five-height local release promoted to a
   first-class candidate type), so renegotiations that matter clear the
   bar as units and stay dated in the ledger.

## Current recommended self-child configuration

`DMESH_SELF_CHILD=1 DMESH_REFRESH_EVERY=1` (all other flags as shipped).
State: -0.011 on the book excluding the issuance-boundary corner, +0.031
overall under the writeup's WALA-0-exclusion protocol, zero null-battery
admissions, no topology-collapse mode, ~4x wall-time, stable dated
certificates. Remaining deficit is localized and mechanistically understood.

## Addendum 2: upward-flow experiments (all MEASURED, seed 7, all +rf1)

Hypothesis (correct, per the corner evidence): the baseline's residual edge
is Gauss-Seidel's upward information flow - child evidence continuously
updates ancestors, keeping every coefficient near E[.|current model]
(martingale property); freezing accumulates the upward correction as
sub-threshold drift. Implementations tested:

| variant | held-out dev | verdict |
|---|---|---|
| baseline joint solve | 6.455 | reference |
| strict freezing (plain rf1) | 6.486 | best self-child |
| reneg-1: parents thawed, prior at no-change | 6.595 | two-rooms damage |
| reneg-1 + grandparent anchor (`DMESH_SC_RENEG_ANCHOR=1`) | 6.497 | anchor recovers 0.10; still no win |
| reneg-2 (one-ring) + anchor | 6.586 | worse |
| upward pass v1: prior-consensus BP, ladder tau (`DMESH_SC_UPWARD=1`) | 8.575 | tau collapse reads child signal as parent error |
| upward pass v2: empirical slab, damped, consensus>=2 | 7.663 | no data term: pass fights the likelihood, never converges |

Two mechanistic lessons, both now measured twice from different directions:
(1) dropping the grandparent anchor from any renegotiation costs ~0.1
(the hierarchical target IS load-bearing information about the parent);
(2) an upward channel without the data term optimizes prior-consistency
against the likelihood - where child deviations encode true surface shape,
the "consensus" never re-centers and the pass drags ancestors away from the
data. The exact version of upward flow (grandparent prior x own data x
child messages) is Gaussian BP whose fixed point is the joint MAP - i.e.,
the baseline. Partial upward channels inject inconsistency at their
boundary that costs more than the drift they repair.

Standing interpretation: full-joint (upward flow everywhere, no decision
identity) and strict freezing (identity, pays ~+0.03 renegotiation debt)
are the two coherent points measured so far. The one untested middle path
with a principled boundary remains the GATED group move: price a family's
joint re-solve as one decision, so a renegotiation happens only when it
pays for its own boundary damage on evidence. Both experimental flags are
retained as measured negatives.

## Addendum 3: joint solve + self-child composition (MEASURED, seed 7)

`DMESH_SC_JOINT` decouples the self-child decision architecture (gated
delta candidates) from the frozen-increment value semantics: candidates and
gating unchanged, heights re-solved jointly every round as in the baseline.

| variant | held-out dev |
|---|---|
| baseline | 6.455 |
| frozen self-child + rf1 | 6.486 |
| joint + self-child, no-change priors (`SC_JOINT=1`) | 6.659 |
| joint + self-child, hierarchical priors (`SC_JOINT=2`) | 6.505 |

Both outcomes were predicted before running and confirm a structural fact:
a self-delta's evidence IS the coefficient's distance from the joint
optimum. Under continuous joint re-solving that distance is consumed before
the gate sees it - with hierarchical prior centers the self-deltas are
near-vacuous and the mode reduces to the baseline plus interference (null
population flooded with near-zero scores; tau contaminated by micro-deltas:
the 0.05 residual); with no-change centers the only remaining "evidence" is
prior tension, and the gate becomes a de-shrinkage pump (1,181 admissions,
worst held-out of the structured variants). The self-child gate and the
joint solve are complementary, not composable: the gate is only meaningful
over values that are NOT continuously re-optimized. Decision-identity
architectures must therefore either freeze values (the original
formulation) or attach identity to something other than coordinate values -
the innovation ledger (per-admission full-field deltas over a jointly
solved state), which remains the designed next step alongside the gated
group move.

## Addendum 4: gated-measurement smoother (MEASURED, seed 7, +rf1)

`DMESH_SC_SMOOTH`: between rounds, rendered heights solve the Gaussian
system {hierarchical prior} + {gated measurements only} (each admission and
base-fit coefficient recorded as a noisy observation of its height). This
is the correctly-posed version of the upward-flow idea: grandparent prior
pulls the parent, child MEASUREMENTS inform ancestors through the prior
coupling, raw data enters only through the gate, and observations are
distinguished from prior residuals (the defect that made the addendum-2
upward passes diverge).

| variant | held-out dev |
|---|---|
| strict freezing (+rf1) | 6.486 |
| smoother v1 (penalized obs, slab from smoothed field) | 6.666 |
| smoother v2 (algebraic un-shrunk obs, slab from measurements) | 6.716 |
| smoother v3 (penalized obs, slab from measurements) | 6.558 |

Isolations: estimating the slab from the smoothed field is a ratchet
(shrink -> smaller slab -> more shrink; fixing it recovers 0.11); the
algebraic un-shrink h + lambda*sigma2*(h - m) is ill-conditioned at weak
coordinates and poisons the measurements (-0.16). Best smoother still
loses to freezing by 0.07.

## Synthesis: why value-side upward flow keeps losing

Nine variants of "let information move frozen values outside the gate" are
now measured (reneg x3, upward-BP x2, joint x2, smoother x3 minus overlap),
every one worse than strict freezing. The pattern has one consistent
reading: all of these channels move values through the PRIOR system, and
the depth-tier exchangeable prior is not accurate enough to be load-bearing
as an information channel - which is the writeup's own EB conceptual audit
conclusion reached from an independent direction (depth is the wrong
exchangeability unit; the learning signal is corrupted; one tau spans
Keff-0.01 corridors and Keff-100 workhorses). The baseline's upward flow
works because it flows through the joint DATA likelihood; the gate forbids
frozen values from touching data outside admissions; so the self-child
architecture's renegotiation debt (~+0.03) is, at current prior quality,
irreducible by value-side mechanisms.

Consequences for the roadmap, in order of leverage: (1) improve the prior
itself - the spike-and-slab per (depth x evidence-geometry) cell from the
writeup's EB endgame - and only then retry the measurement smoother, whose
machinery is correct and in the tree; (2) gated group moves (data-driven
renegotiation that pays FDR); (3) the innovation ledger (joint values,
decision identity attached to per-admission full-field innovations rather
than coordinates).
