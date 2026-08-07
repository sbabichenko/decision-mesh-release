# Reproducibility guide

## Contents

- `experiments/controlled_null_decomposition/frozen_inputs/`: exact controlled-null CSVs used in the report.
- `experiments/controlled_null_decomposition/generate_controlled_nulls.py`: documented mechanism generator for new deterministic draws.
- `experiments/controlled_null_decomposition/run_report_battery.sh`: exact report matrix.
- `experiments/controlled_null_decomposition/run_one.sh`: one-run environment and command.
- `experiments/controlled_null_decomposition/analyze_one.py`: metric extraction.
- `experiments/corrected_homogeneous_certified_null/`: exact five HCN inputs, generator, and before/after summaries.
- `manifests/frozen_run_manifest.csv`: one row per frozen completed run, including input, mode, stopping rule, parent scale, and output metrics.
- `manifests/scenario_manifest.csv`: exact input hashes and realized rates.

## Build and test the repaired code

```bash
cmake -S source -B build -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The production source directory in this delivery is intentionally only the files touched by the repair plus its CMake file; apply `lindsey_guard.patch` to the full reviewed tree for a complete build context.

## Run the controlled report battery

```bash
experiments/controlled_null_decomposition/run_report_battery.sh   /path/to/decision_mesh   experiments/controlled_null_decomposition/frozen_inputs   /tmp/null_decomposition_runs
```

This is computationally nontrivial. Individual rows can be rerun with `run_one.sh`; the exact arguments are in `manifests/frozen_run_manifest.csv`.

## Regenerate mechanism-equivalent controlled scenarios

```bash
python experiments/controlled_null_decomposition/generate_controlled_nulls.py   --design experiments/controlled_null_decomposition/instruments/ginnie_design.csv   --truth-beta experiments/controlled_null_decomposition/instruments/null_truth_beta.npy   --out /tmp/new_null_draws --seeds 0 1 2
```

These are fresh deterministic draws of the declared mechanisms. Use the packaged `frozen_inputs` to reproduce the published numerical tables exactly.

## Reproduce HCN inputs exactly

```bash
python experiments/corrected_homogeneous_certified_null/generate_hcn_corrected_portable.py   --out /tmp/hcn --seeds 0 1 2 3 4
```

The generated HCN files should match the packaged frozen inputs byte-for-byte under the listed NumPy/Pandas serialization behavior; `manifests/scenario_manifest.csv` contains the acceptance hashes.
