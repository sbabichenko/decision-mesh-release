# Build validation - July 30, 2026

A clean Release build was configured with GCC 14.2.0 using the packaged
`scripts/build_release.sh` script.

Results:

- `decision_mesh`: built successfully;
- `test_lindsey_guard`: built successfully;
- CTest: 1/1 tests passed;
- compiler warnings promoted by the existing build configuration: none fatal.

The Ginnie NLL comparison and full-objective fixed-point audit are preserved in
`docs/coupled_penalty_validation/`. Those artifacts were produced on the benchmark
input used for the 2.9410 / 2.9520 / 2.9443 comparison.
