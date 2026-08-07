# Phase B1 implementation status

Status: **complete** (2026-07-28).

Implemented and regression-tested:

- one-dimensional exact penalized candidate profiles;
- shared inactive-candidate and active hierarchical-column construction;
- live re-profile immediately before commitment;
- negative-gain admission exclusion;
- full-objective commit invariant;
- dense one-dimensional profile validation;
- shared active hierarchical design engine;
- hierarchical support diagnostics before EB scale estimation;
- deepest-first fixed-topology retirement on raw training columns;
- per-batch zero-training-change assertion before re-optimization;
- thin/empty candidate classification;
- named regime-straddle logging and grid check;
- historical stale-vertex census audit;
- deterministic clean-run comparison;
- eruption regressions for seeds 202 and 404.

Headline closure numbers:

- original held-out NLL: 2.9335;
- original WALA 0-2 fit: -6.96 against law-marginal -7.95;
- exact candidate gains: 0 negative, 24/24 grid checks passed;
- commit checks: 4,029 passed, 0 failed;
- stale census: maximum standardized residual move 1.07e-06;
- seed 202 NLL: 2.9498, no eruption;
- seed 404 NLL: 2.9450, no eruption;
- deterministic model artifacts: byte-identical.

See `B1_TEST_RESULTS.md` for the full closure record.

Next phase:

- B2 chooses the exact-gain gate statistic;
- B3 recalibrates its null;
- B4 runs v3/v4 certification acceptance;
- B5 runs the complete four-part acceptance battery;
- B0 scaffolding removal is repeated after the legacy gate path is gone.
