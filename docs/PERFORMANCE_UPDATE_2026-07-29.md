# Decision Mesh performance release - 2026-07-29

## Frozen benchmark

The benchmark uses the included `data/ginnie_design.csv` (120,595 Ginnie design rows), an even/odd held-out split, `q=0.10`, 24 points per face, exact likelihood fitting, B3 signed-root scoring after the stage-0 pilot, and `DMESH_GATE_EXTRA_VAR=0.35`. The exact command is encoded in `scripts/run_ginnie_benchmark.sh`.

On the release machine, the clean Release build produced:

| Gate | External wall time | Program `TOTAL` | Final active vertices | Held-out deviance/pool | X2/info |
|---|---:|---:|---:|---:|---:|
| Direct candidate BH | 7.43 s | 7.1 s | 1,261 before final retirement | 6.687 | 0.1493 |
| Adaptive Lindsey | 6.51 s | 6.2 s | 1,055 before final retirement | 6.683 | 0.1510 |

These are full-pipeline times, including final retirement/refitting, leverage correction, PM diagnostics, and CSV export. The stage-2 adaptation timestamps were about 4.9 s for BH and 4.5 s for Lindsey. Earlier conversational references to a “4/7 build” mixed the adaptation timestamp with full executable wall time; this release records both scopes explicitly.

## Accepted implementation changes

- Lindsey Poisson fits fail closed if IRLS does not converge or the fitted density fails mass, finiteness, or concavity checks.
- The pipeline stops exactly when the stage variance input and output are unchanged, avoiding a deterministic repeated third stage.
- Release builds use `-O3`, optional native tuning, and LTO when supported.
- The exact height solver materializes sparse likelihood rows once per solve and reuses one visit buffer across coordinate updates.
- Counts, exposures, frozen likelihood tempers, and sparse column values are stored in a canonical solver-row representation.

## Reproducibility

Run:

```bash
./scripts/build_release.sh
./scripts/run_ginnie_benchmark.sh ./build-release/decision_mesh bh
./scripts/run_ginnie_benchmark.sh ./build-release/decision_mesh lindsey
```

The release includes reference mesh and hierarchy files. `scripts/verify_reference.sh` rebuilds both variants and checks their SHA-256 hashes.

## Rejected experiments

The release excludes optimizations that were fast but changed adaptive results, including persistent candidate-column caches with incomplete invalidation and descendant-only column traversals that did not reproduce the canonical targeted columns. It also excludes performance-neutral or slower workset, heap-frontier, warm-start, and scratch-workspace experiments.
