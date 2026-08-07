# Decision Mesh

Deterministic adaptive piecewise-linear regression on a conforming triangular mesh, with empirical-Bayes shrinkage, an empirical local-FDR refinement gate, and residual-variance diagnostics.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Run the simulation pipeline:

```bash
./build/decision_mesh [sigma_u=0.30] [q=0.10] [scale=1] [points_per_face=12] [seed=7]
```

The implementation is serial. Stable geometric ordering and fixed probe sequences preserve deterministic model outputs.

## Source layout

- `gated_main.cpp`: argument parsing, dataset construction, and top-level error handling only.
- `pipeline.{h,cpp}`: deterministic fitting orchestration and refinement stages.
- `run_config.{h,cpp}`: grouped simulation, fitting, diagnostic, and output configuration.
- `dataset.{h,cpp}`: simulation and checked CSV loading.
- `output.{h,cpp}`: checked CSV exports using `std::ofstream`.
- `diagnostics.{h,cpp}`: leverage correction, sparse CG, and Paule--Mandel diagnostics.
- `lfdr.{h,cpp}`: Lindsey density estimation and local-FDR calculations.
- `mesh.*`, `vertex.*`, `edge.*`, `face.*`, `tree.*`: adaptive mesh implementation.
- `units.h`: explicit conversions between natural log-odds and the historical centi-log-odds storage scale.

## Common environment options

| Option | Purpose |
|---|---|
| `DMESH_DATA=<csv>` | Load real data with columns `wala,wac,n0,k_term`. |
| `DMESH_CHIRP=1` | Use the multi-scale chirp target. |
| `DMESH_AMP=<value>` | Set synthetic-surface amplitude. |
| `DMESH_PERM=<0..1>` | Persistent share of simulated pool variance. |
| `DMESH_ATTR=<value>` | Simulated attrition coefficient. |
| `DMESH_NOISESEED=<seed>` | Resample outcomes while holding locations fixed. |
| `DMESH_NOEB=1` | Disable empirical-Bayes shrinkage for ablation. |
| `DMESH_SPLIT_MODE=local` | Default real-data holdout: keep pools intact and balance them within local WALA/WAC blocks. |
| `DMESH_SPLIT_MODE=parity` | Reproduce the historical even-row/odd-row split for comparison only. |
| `DMESH_SPLIT_BLOCK_SIZE=<4..256>` | Maximum number of pools per local split block; default 24. |
| `DMESH_METRIC=1` | Enable the learned separable axis map. |
| `DMESH_ORACLE=<vertices>` | Run the noiseless oracle placement experiment. |
| `DMESH_DUMP=<prefix>` | Write surface, observation, mesh, and pool CSV outputs. |
| `DMESH_LEVERAGE_PROBES=<1..64>` | Set the deterministic leverage-probe count; default 6. |
| `DMESH_TAU0_LOGIT_SD=<sd>` | Override the weak default EB prior SD. |
| `DMESH_TAU_FLOOR_LOGIT_SD=<sd>` | Override the EB prior-SD floor. |
| `DMESH_TRUNC=1` | Opt into the experimental truncation-corrected tau moment. |

Diagnostic-only switches remain in `run_config.h`; they are intentionally centralized rather than read ad hoc throughout the pipeline.

The split membership is exported as `is_train` in `run_obs.csv`. The local split uses geometry and exposure only; outcomes are not used. See `LOCAL_BALANCED_SPLIT_BOOTSTRAP_2026-07-30.md` for the same-file parity comparison and the 20-seed audit.

## Removed compatibility paths

The parent-propagating EB prior, all-active tau estimator, and repeated-trace PM solver were removed. Their former switches (`DMESH_LEGACY_EB`, `DMESH_TAU_ALL_ACTIVE`, and `DMESH_PM_LEGACY`) now produce an explicit error rather than silently selecting stale code.


## Current synthetic benchmark (July 2026)

The default synthetic generator now uses a controlled rectangular design rather than the earlier cohort wedge:

- no attrition;
- every simulated pool is observed for WALA 1--60;
- WAC is sampled independently over its full range;
- displayed horizontal coordinate is `x = WALA / 12`;
- the truth is a two-scale exponential-cosine surface whose horizontal phase is warped by
  `z_warp = (0.2 + 0.8*z)*z`, so frequency increases toward larger WALA.

This controlled benchmark is intended to test geometric adaptation without confounding it with support artifacts. Real data continue to be loaded with `DMESH_DATA=<csv>` and are normalized to the observed WALA/WAC rectangle.

## Real-data hierarchy diagnostics

`DMESH_DUMP=<prefix>` exports the final mesh and observation files. The current `output.cpp` additionally exports hierarchy metadata needed to inspect deep refinement chains and empirical-Bayes behavior. The July 2026 Ginnie Mae experiments found that very fine clusters can be undershrunk: manually smoothing the finest local mesh regions improved held-out binomial deviance, while a first broad local-branch EB model overshrank healthy regions and is therefore not part of the recommended estimator.

## Geometry/coefficient lifecycle and boundary re-audit (July 2026)

The current branch separates mesh conformity from statistical freedom.

- An NVB completion midpoint is created geometrically with zero surplus, `height = mu_lin()`, and does not become a free coefficient merely because conformity requires it.
- A geometry-only midpoint remains eligible for later gate testing and may be promoted without changing the fitted surface at the instant of promotion.
- A previously admitted non-root coefficient with zero current training curvature becomes dormant during adaptation, keeps its last accepted height, and wakes if support reappears.
- After topology is fixed, still-zero-support coefficients are retired to zero surplus and the remaining supported coordinates are swept without any later refinement.
- Free nodal coordinates with fewer than 10 current training observations or effective support below `DMESH_MIN_CURRENT_KEFF` (default `0.5`) are held at their last accepted values as a temporary containment rule.
- `DMESH_HEIGHT_LIMIT_LOGIT=<L>` enables a diagnostic absolute height guard of `L` logits. The default is no height limit. With the lifecycle and support fixes, matched 12-logit, 15-logit, and unlimited real-data runs converge to the same practical fit.

The former hard-coded five-logit ceiling was removed. On the shipped real-data snapshot it mechanically imposed a minimum predictor near -3.704 logits and accounted for most of the previously reported WALA-zero under-dive.

## Phase B1 hierarchical support semantics

`hierarchical_design.{h,cpp}` is the single active-column builder for the exact
stage solver and final support cleanup. Final retirement is based on the raw
hierarchical column on explicit training rows, proceeds deepest-first, and
asserts zero immediate training-prediction movement. Nodal-star emptiness no
longer retires or dormants a coefficient that remains supported through
constrained descendants. See `DORMANT_RETIREMENT_CHANGE.md` and
`B1_IMPLEMENTATION_STATUS.md`.

## Phase B0 post-B1 cleanup

After exact live candidate profiling and commitment were completed in B1, the
legacy saturation-triggered exact fallback was regression-inert on the original
split, eruption seeds 202/404, and the 20-split stability battery. It and its
admission-pause state were removed. The 12-logit height guard and four-logit
IRLS working-step cap remain required by the still-legacy gate path and are
retested after gate-statistic migration. See `B0_TEST_RESULTS.md`.

## Phase B2 exact signed-root gain

Every candidate now exposes the direction-preserving exact penalized-gain score

`sign(delta_hat) * sqrt(2 * [Q(delta_hat) - Q(0)])`.

Use `DMESH_B2_SCORE_DUMP=<csv>` to write the shadow ledger and
`b2_analyze.py <scores.csv> --out <dir>` to summarize objective identities,
legacy/exact ranking disagreement, probability regimes, and roundwise scale.
The production gate is unchanged.  `DMESH_B2_USE_SCORE=1` is an untuned
research diagnostic only; the complete-null experiment shows that the raw
score is stage dependent and must be recalibrated before certification claims
can transfer.  See `B2_TEST_RESULTS.md`.
