# Decision Mesh writeup revision - posterior pool effects

Base source: `DecisionMesh_writeup_current_benchmark_reruns_2026-07-31.tex`, recovered from the user's prior file library.

## Changes

- Updated the document date to August 1, 2026.
- Added the explicit three-way decomposition:
  existing-model offset + shared WALA/WAC surface + posterior pool effect + remaining noise.
- Revised the single-snapshot identification language. A single cross-section does not identify whether heterogeneity is temporary or persistent, but it does support one-shot empirical-Bayes posterior pool effects.
- Added a real-data subsection on posterior pool heterogeneity for all 120,595 Ginnie Mae pools.
- Added four corrected figures:
  - posterior pool effects across the full WALA-WAC support;
  - all-pool histogram with a symmetric-log x-axis;
  - central 99.8% histogram with an approximate percentage-change axis;
  - exposure-weighted central histogram.
- Corrected the plotting units: the C++ `run_pools.csv` output stores `u_hat` in centi-log-odds, so all displayed log-odds effects divide `u_hat` by 100.
- Added an accurate summary of exploratory shrinkage results: direct surplus-neighbor pooling, WAC-strip variance sharing, and full graph-hierarchy replacement were not retained; graph patches were modestly positive; crossed function-level covers were strongly positive but remain experimental.
- Strengthened the discussion and conclusion around the decomposition of shared geometry, pool-specific posterior heterogeneity, and remaining noise.
- Left the shipped code release unchanged.

## Figure conventions

- The scatter plots all 120,595 pools. Only the color scale is clipped at the 1st and 99th percentiles.
- The full histogram includes every posterior effect.
- The central histogram displays the central 99.8% for readability and states that the full histogram includes all pools.
- The top percentage axis uses the first-order approximation `100 * u` percent and is described as an approximate change in termination odds for modest effects.
