# Experimental pooled-permutation admission gate

Status: implemented and tested, but **not adopted as the release default**.

The production source remains in `../../core/`. This directory contains a full
experimental source tree plus a patch against the production source.

## What is implemented

- candidate-specific permutation moments from the actual hierarchical star;
- 32 deterministically seeded permutation draws per candidate by default;
- cross-star pooling of standardized reference draws by star-size band;
- direct BH (`DMESH_PERM_BH=1`) or theoretical-null lfdr;
- configurable effective-support threshold (`DMESH_PERM_MIN_KEFF`);
- a second effective-support check at candidate commit;
- unit tests for the permutation-moment implementation.

This is not yet the memo's deterministic double-saddlepoint tail calculation.

## Build and test

```bash
cmake -S source -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The packaged clean build passed both tests. Logs are in `cmake.log`,
`build.log`, and `ctest.log`.

## Experimental operating point

The best heavy-tailed certified-null setting was:

```bash
export DMESH_PERM_GATE=1
export DMESH_PERM_BH=1
export DMESH_PERM_DRAWS=32
export DMESH_PERM_MIN_KEFF=4
# pass q=0.15 to the executable
```

It made no certified-window admissions in 20 heavy-tailed seeds and improved
outside-window fit, but it underfit the real Ginnie cross-section. After
removing WALA-zero pools and refitting ten matched splits, mean held-out
deviance was 6.8108 versus 6.4968 for the shipped gate. Therefore the branch
is retained for research, not silently enabled.

## Evidence

See `../../results/permutation_gate/` for:

- the 34-run complete-null battery;
- the Gaussian corrected certified null;
- the 20-seed heavy-tailed certified null;
- the unfiltered real-Ginnie comparison;
- the real-Ginnie refit after removing WALA zero.

The original independent memo is preserved as
`PROPOSAL_MEMO_FROM_OTHER_CHAT.md`.
