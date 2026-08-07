# Decision Mesh complete release - updated 2026-08-01

This package combines the stable C++ release, the adopted post-release
quadratic inherited-child update, current benchmark evidence, the real Ginnie
Mae posterior-pool-effect run, and the August 1 permutation-gate audit.

The pooled-permutation gate is included as a separate experimental source tree.
It is not enabled in `core/`: it passed the complete-null and heavy-tailed
certified-window batteries, but underfit the real Ginnie cross-section.

## Adopted estimator in this package

- exact coupled hierarchical penalty solver;
- locally balanced train/holdout pool split;
- hierarchical MAP fit;
- quadratic inherited-child prior target (`DMESH_QUAD_PRIOR=1`);
- Lindsey empirical-null gate when the `lindsey` benchmark mode is selected.

The quadratic update was applied to the 2026-07-30 local-split release using
`experiments/post_release/quad_prior_coupled.patch`. The geometry-patch,
crossed-cover, two-scale, and interaction-precision files are retained under
`experiments/post_release/` as experiments; they are not silently enabled in
the current release path.

## Build

```bash
cmake -S core -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release -j
ctest --test-dir build-release --output-on-failure
```

## Run the adopted current configuration

```bash
./scripts/run_current_model.sh ./build-release/decision_mesh 7 output/current_seed7
```

For the bundled Ginnie design and Lindsey empirical-null gate:

```bash
./scripts/run_ginnie_benchmark.sh ./build-release/decision_mesh lindsey output/ginnie_lindsey 7
```

The Ginnie script enables the quadratic inherited-child target.

## Important directories

- `core/`: C++ estimator and Python analysis utilities.
- `scripts/`: build and benchmark scripts.
- `data/`: bundled Ginnie design and test data.
- `writeup/permutation_gate_audit/`: latest PDF, LaTeX, figures, and summary tables.
- `writeup/posterior_pool_effects/`: previous August 1 writeup revision.
- `results/current_benchmarks/`: July 31 controlled rerun summaries.
- `results/ginnie_posterior_pool_effects/`: real-data run outputs and corrected figures.
- `experiments/pooled_permutation_gate/`: full experimental C++ branch, patch, tests, build logs, and proposal memo.
- `results/permutation_gate/`: complete-null, certified-null, heavy-tail, real-data, and WALA-zero-exclusion results.
- `experiments/post_release/`: other promising or negative experimental branches not adopted by default.

## Pool-effect units

`run_pools.csv` stores `u_hat` in centi-log-odds. Convert to log-odds with
`u_hat / 100`. The corrected figures in this release use that conversion. The
secondary percentage axis is the small-effect approximation `100 * u`,
interpreted as an approximate percentage change in the constructed attrition-proxy odds.


## Ginnie outcome qualification

The bundled April `ginnie_design.csv` does not contain SMM. Its response was
constructed as `round(Original Aggregate / AOLS) - current loan count`. Because
AOLS describes the currently disclosed pool composition, this ratio is not an
exact original loan count. The latest writeup therefore calls the response a
constructed loan-count attrition proxy, removes WALA-zero pools in the main
sensitivity audit, and does not interpret those probabilities as SMM.
