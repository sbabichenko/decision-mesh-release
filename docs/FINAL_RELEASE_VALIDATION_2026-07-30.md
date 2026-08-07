# Final local-split release validation — July 30, 2026

The final archive was extracted into a clean directory and built with:

```bash
./scripts/build_release.sh /tmp/decisionmesh-build
```

The `lindsey_guard` CTest passed (1/1). The extracted release binary was then run on the shipped `data/ginnie_design.csv` with the default local split and seed 7:

```bash
./scripts/run_ginnie_benchmark.sh /tmp/decisionmesh-build/decision_mesh bh /tmp/seed7 7
python3 core/eval_harness.py /tmp/seed7/run_obs.csv data/ginnie_design.csv core/auditing/law_shared.json
```

Reproduced results:

- marginal NLL: 2.7129;
- ramp / mid / seasoned NLL: 2.665 / 3.046 / 2.376;
- split: 60,250 training rows and 60,345 held-out rows;
- final topology: 967 faces and 526 active vertices;
- internal executable runtime: 6.4 seconds on the validation host.

The 20-seed audit in `local_split_bootstrap/` uses seeds 101--120. Its primary mean includes the two topology-collapse runs (seeds 101 and 107).
