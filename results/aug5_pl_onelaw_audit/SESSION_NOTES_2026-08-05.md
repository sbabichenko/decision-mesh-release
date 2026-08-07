# August 5, 2026 session - exact-solver default, one-law gate, converged law

Full account in writeup/permutation_gate_audit/DecisionMesh_writeup_aug5_audit_2026-08-05.{tex,pdf},
Section "August 5 audit". Summary:

- External diagnostic's eruptions reproduced and traced to DMESH_PL_SOLVER
  being off in its harness. With the exact solver on, eruptions vanish on all
  ten MC masks and the mesh beats the tensor GAM 10/10 (dev 1.175 vs 1.243).
- Pooled-permutation gate retested on SMM under the exact solver: null control
  and fit both fail in this regime; rejection stands and strengthens.
- One-law gate activated; joint surface-law fixed point converges (no
  oscillation). Converged SMM pool law per stratum (pi_heavy, tail var, total;
  excess kurt): (3.7%,1.36,0.143;8.5) (7.6%,1.62,0.242;8.1)
  (1.8%,0.31,0.029;5.4) (1.0%,0.004,0.004;0.0). The oscillation-era law is
  retracted (totals inflated 4-17x by surface-error absorption).
- Two null-generator artifacts diagnosed and fixed (stratified Jensen shift;
  central-matching vs composite pricing -> empirical-null sd floor).
  Corrected flat-marginal battery: one-law+floor 0.2-0.3 false adm/run,
  naked Lindsey 11/run with fat seed tail.
- DMESH_LAW_WEIGHTS convicted by ablation (double-charges pool noise with the
  one-law gate; +0.16 dev on the regressing split); removed from the
  recommended config. DMESH_LAW_WEIGHTS_FINAL (two-currency GLS final solve)
  implemented, measured (-0.0004 mixture NLL, +0.004 deviance), NOT adopted:
  complexity not justified. Experimental flag only.
- Moment instruments calibrated: clean-null falsification exact; PM/energy
  run ~3.5x hot under heavy tails; EM mixture law is the primary variance
  instrument; pairwise sigma_a^2 needs the two-month panel.

RECOMMENDED CONFIG:
  DMESH_PL_SOLVER=1
  DMESH_ONELAW="0.1428,0.2415,0.0286,0.0038;8.5,8.1,5.4,0.0"
  DMESH_SIGMA0_MIN=1.0
Measured: MC-10 held-out deviance 1.189 (GAM 1.243; naked Lindsey 1.175),
mixture NLL 0.930, ~148 admissions/split, 0.2-0.3 false admissions/run on the
honest null battery.

CODE CHANGES IN THIS TREE (vs 2026-08-01):
  core/dataset.cpp  - DMESH_SPLIT_MASK external mask replay (also applied in
                      experiments/pooled_permutation_gate/source/dataset.cpp)
  core/pipeline.cpp - DMESH_LAW_WEIGHTS_FINAL experimental flag (frozen-
                      topology GLS temper; not part of the recommended config)
SHA256SUMS.txt predates these changes and is stale for the two files above.

DATA/ARTIFACTS HERE:
  ../..//data/smm_design_202606.csv       - the June SMM design (33,740 pools)
  mc10_split_masks/                       - the ten external MC masks (rep_06_replacement = seed 1111)
  *.csv                                   - all result tables from the session
  figures/                                - annotated star plots (shipped vs PL, human-scale, pool-effect coloring, high-Keff, worst-pool, sparse region)
  decision_mesh_star_inspector.html       - interactive split-1 star dashboard (both arms embedded)
OPEN ITEMS: admission stability / cross-fitted admission layer (rep-8 class),
gate_k4 CF term wiring, cross-split law confirmation, two-month panel.

## Evening addendum - ridge counterfactual and marginal height solve

- Ridge counterfactual (v1894 star, split 1): held-out mixture ll peaks at
  ~40% of the fitted surplus (peak ~8%, not 23.4%); the binomial height solve
  cannot trade surface for pool effect (u==0 in its objective).
- marginal_final.py: frozen-topology law-marginal final pass (envelope-theorem
  score at the posterior mode, Laplace-tempered curvature, early stop at the
  marginal-ll peak). Ridge lands ~10%; blue one-sided cluster dissolves.
  'all' mode = post-hoc all-data oracle at fixed mesh.
- Balanced-split check (DMESH_SPLIT_MODE=local, seeds 1-10): certification is
  nearly free (rec vs naked dNLL +0.0004); marginal pass NLL 0.92423->0.92201
  (better 7/10). External masks kept for diagnostic comparability.
- Writeup updated in place: star-diagnostics subsection, ridge counterfactual,
  marginal solve, oracle figures. PDF now 89pp, compiles clean.
- New results here: marginal_solve_mc10.csv, balanced_splits_nll.csv,
  marginal_final.py, figures/ additions (offending star x3, 2x2, honest surface).

## Temporal-layer design (added to writeup, sec "Outlook: the temporal layer")

Meshes stay adaptive per month; months couple through the weak formulation
(test-function moments of each P1 surface, exact on any mesh). Dynamics = a
simple linear second-order parabolic PDE with O(10) operator coefficients,
GLS-estimated with weights from the fits' posterior covariances. Advection
b=(1,0) is KNOWN (WALA ages deterministically) -> null model is pure
transport (surface constant in vintage coordinates); issuance is the inflow
boundary. Test-dictionary choice = the temporal admission problem (FDR over
test functions); dictionary must be pre-registered before any temporal fits.
Status: design only, nothing measured.

## Admission churn decomposition (late night)

Sam measured churn on the ten balanced seeds (determinism confirmed: same mask,
different seed = bit-identical). Decomposition against two gate-correct
simulated worlds (full-data rec surface + converged law): real evidence-wtd
Jaccard 0.386 vs sim 0.51 -> ~79% of churn is the power floor; path-dependence
share of absences 20% vs 8-9%; surface tail mass (sd>0.3) 4.9% vs 0.8%.
Remedy benchmark is the sim floor, not 1.0. Two independent real-seed
measurements disagree on level (J 0.27 vs 0.53) - metric reconciliation
pending. Data: churn_decomposition.csv; sim worlds simworld{1,2}.csv from
full_rec surface. Writeup: new subsection sec:aug5-churn + Discussion item 2.
