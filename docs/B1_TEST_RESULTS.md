# Phase B1 closure results

Date: 2026-07-28

## Configuration

```text
DMESH_PL_SOLVER=1
DMESH_HIER_FIT=2
DMESH_HEIGHT_LIMIT_LOGIT=12
DMESH_IRLS_STEP_CAP=4
DMESH_HIER_MAP=1
DMESH_SPLIT=1
DMESH_B1_VALIDATE=1
DMESH_B1_COMMIT_VALIDATE=1
DMESH_B1_STALE_VALIDATE=1
```

The height limit and IRLS step cap remain enabled because the legacy gate path
still uses them. The saturation-triggered exact fallback was removed in the
post-B1 B0 retest: exact live candidate re-profiling made it regression-inert.

## Exact candidate profile and commit invariants

On the original split:

- 24,461 candidate profiles attempted;
- 24,412 scored;
- 49 excluded because the training hierarchical column was empty;
- 0 tempered-out nonempty columns;
- 0 negative exact gains;
- 24/24 dense one-dimensional grid profiles passed;
- 4,029 accepted commits checked against the live full penalized objective;
- 0 commit-check failures;
- maximum stored-versus-realized gain error: 2.33e-10 nats;
- maximum stored-versus-realized optimum error: 1.87e-08 centi-logits.

The inactive-candidate and active-coordinate column builders are shared. A live
re-profile is performed immediately before commitment so the committed
coordinate is priced in the same objective and hierarchy that it receives after
activation.

## Regime-straddle regression

The nearest actual admission candidate to raw `(WALA,WAC)=(0,4.52)` is:

```text
normalized location: (0.015625, 0.390625)
null:               -3.805815 logits
exact optimum:      -5.574881 logits
exact gain:          4.82920046 nats
```

The named fixed-topology pin at normalized `(0,0.375)` also passed its dense
grid profile check.

## Hierarchical retirement

On the original split:

- 19 genuinely empty hierarchical columns retired in four deepest-first batches;
- 14 nodally empty but hierarchically supported coefficients woken;
- maximum immediate training-prediction change from every retirement batch: 0;
- maximum immediate held-out change: 0.963 centi-logits;
- no free coefficient is retired or dormant after cleanup.

The later fixed-topology EB re-estimation and exact solve is reported separately
from retirement. It is allowed to move predictions because it solves a new
penalized objective after empty coordinates have been removed.

## Historical stale-vertex census

The saved census contains 1,369 historical rows. The current topology matches
1,299 by exact normalized coordinates; 831 matched rows remain free and are
directly auditable. Of those, 237 had historical `|z_gap|>2`.

At the current exact coordinate fixed point:

- maximum absolute coordinate move: 4.42e-05 centi-logits;
- RMS coordinate move: 1.56e-06 centi-logits;
- maximum absolute standardized move: 1.07e-06;
- rows with current `|z_move|>0.01`: 0;
- maximum residual exact gain: 1.75e-10 nats.

## Determinism

Two independent clean original-split runs are byte-identical for:

- observation predictions;
- topology;
- hierarchical vertex dump;
- surface dump;
- pool dump;
- exact-gain dump;
- stale-census dump.

Only elapsed-time text differs in stdout/stderr.

## Predictive regressions

| run | held-out NLL | law-marginal WALA 0-2 dot | fitted WALA 0-2 | eruption |
|---|---:|---:|---:|---:|
| original split | 2.9335 | -7.95 | -6.96 | no |
| seed 202 | 2.9498 | -8.33 | -6.99 | no |
| seed 404 | 2.9450 | -8.10 | -6.72 | no |

The original split satisfies both Phase-B predictive pins: NLL is within 0.01
of 2.9347 and the fitted WALA 0-2 base is at most -6.9. Seeds 202 and 404 are
regression tests for the former eruption class, not separate boundary-acceptance
cells.

## B1 decision

B1 is complete. Exact candidate profiling, live commitment, hierarchical
retirement, the stale-vertex regression, the named regime-straddle profile, and
deterministic repeatability all pass.

Not completed in B1:

- migration from the legacy Wald/local-FDR statistic to an exact-gain statistic;
- v3/v4 null recalibration;
- the full 20-split acceptance battery;
- post-migration repetition of B0 scaffolding removal.

## Post-B1 B0 cleanup

The scaffolding-removal test was repeated after exact live candidate commitment.
The saturation-triggered exact fallback was regression-inert and removed; its
20-split ablation produced 0/20 eruptions, mean NLL 2.9449, standard deviation
0.0054, and worst split 2.9537. The height guard and IRLS step cap remain active
in the legacy working-scale gate and were not removed. See `B0_TEST_RESULTS.md`.
