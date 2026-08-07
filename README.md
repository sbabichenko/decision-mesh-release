# August 1, 2026 update

The stable production source in `core/` is unchanged by the pooled-permutation
audit. The latest writeup is in `writeup/permutation_gate_audit/`; the
experimental implementation and its tests are in
`experiments/pooled_permutation_gate/`; and all new batteries are in
`results/permutation_gate/`. The experimental gate is not the release default
because it loses held-out fit on the real cross-section.

The April Ginnie response is now explicitly documented as a constructed
loan-count attrition proxy, not SMM or an exact cumulative-termination count.

# Decision Mesh coupled-penalty + local-split release - July 30, 2026

This archive is the self-contained exact-solver/B3 development release with the July 30 coupled hierarchical-penalty correction and the ordered-holdout repair. The stage solver prices the full Gaussian prior induced by all free descendant surpluses affected by a canonical height move. Real-data holdouts now keep whole pools together inside locally balanced WALA/WAC blocks rather than using CSV row parity.

## Contents

- `core/` - C++17 source, CMake build, tests, and analysis utilities.
- `scripts/` - clean build, benchmark, and reference-verification scripts.
- `data/ginnie_design.csv` - 120,595-row frozen benchmark design.
- `reference/` - byte-pinned BH and Lindsey mesh/hierarchy outputs plus logs.
- `experiments/` - histogram-null and controlled true-null experiment code and compact results.
- `docs/` - implementation handoffs, performance reports, and coupled-penalty validation artifacts.
- `writeup/` - updated LaTeX source, figures, revision note, and compiled PDF.

## Build

Requirements: CMake 3.14+, a C++17 compiler, and standard Unix build tools.

```bash
./scripts/build_release.sh
```

The script configures a Release build, compiles `decision_mesh`, and runs the Lindsey fail-closed regression test.

## Reproduce the frozen benchmark

```bash
./scripts/run_ginnie_benchmark.sh ./build-release/decision_mesh bh
./scripts/run_ginnie_benchmark.sh ./build-release/decision_mesh lindsey
```

The exact environment is encoded in the script, including `DMESH_GATE_EXTRA_VAR=0.35`. The real-data default is `DMESH_SPLIT_MODE=local`; set `DMESH_SPLIT_MODE=parity` only to reproduce the historical ordered-row split. See `docs/LOCAL_BALANCED_SPLIT_BOOTSTRAP_2026-07-30.md` for the same-file comparison and 20-seed audit. Clean-build reproduction is recorded in `docs/FINAL_RELEASE_VALIDATION_2026-07-30.md`.

The `reference/` directory preserves the July 29 pre-correction outputs for historical comparison. Because the coupled-prior correction intentionally changes the fitted fixed point and adaptive topology, `scripts/verify_reference.sh` is a legacy comparator and is not expected to pass against the updated executable.

## Coupled hierarchical-penalty correction

For free-surplus vector `delta = A h`, the prior is

```text
P(h) = 0.5 * (A h - m)' Lambda (A h - m).
```

The exact stage solver now uses the full gradient `A' Lambda (A h - m)` and incrementally updates every affected prior residual during Gauss-Seidel moves. Set `DMESH_COUPLED_FIXED_POINT_AUDIT=1` to print full-objective gradient and boundary-KKT diagnostics after each stage fit.

The historical coupled-penalty comparison in `docs/coupled_penalty_validation/` reports 2.9410 / 2.9520 / 2.9443 for the release baseline, diagonal-only exact stage, and coupled exact stage. Its run log came from an older Ginnie snapshot with pooled rate 0.7851; the frozen CSV shipped here has pooled rate 0.783241. The older input is unavailable, so those numbers are retained as historical validation, not as a reproducible scoreboard for this archive.

On the shipped CSV, seed 7 gives marginal NLL 2.7413 under historical parity and 2.7129 under the locally balanced split. Across seeds 101-120, local-split NLL averages 2.7191 (SD 0.0093); 18 runs cluster at 2.7164 +/- 0.0023, while two topology-collapse runs remain visible and included. Quadratic-prior and ghost-boundary modes retain their specialized treatment and are not covered by the sparse derivative construction.

## Updated paper

- `writeup/writeup_local_split_bootstrap_2026-07-30.pdf`
- `writeup/writeup_local_split_bootstrap_2026-07-30.tex`

The writeup derives the coupled gradient, reports the fixed-point audit, replaces the ordered parity holdout, and records the 20-seed split-sensitivity audit including the two topology-collapse outliers.
