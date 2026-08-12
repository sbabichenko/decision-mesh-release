# Shrinkage prototype: beyond pure depth-based tau

Offline study, on shipped vertex dumps, of moving the empirical-Bayes slab
variance `tau` away from the current *per-depth* exchangeability unit toward
*(depth x evidence-geometry)* cells. Prototyped here first (no C++ loop per
iteration) before any change to `core/mesh.cpp::recompute_tau_sq`.

Prompted by the Aug 6 EB audit's defect (1): one per-depth `tau` spans
K_eff 0.01 corridors and K_eff 100 workhorses. This is the on-ramp after the
`DMESH_FINAL_RESHRINK` TAU_FLOOR bug fix (commit e3bd6f5), so the baseline the
new rule is measured against is the corrected, floor-respecting one.

## Files
- `eb_cells_prototype.py` — recovers raw admission surpluses from a
  `run_hier_vertices.csv` dump, fits `tau` under per-depth / cell / pooled
  schemes, and scores them by held-out-vertices EB marginal likelihood.
- `robust.py` — validation + robustness sweep across seeds (fixed bins,
  tertiles, partial pooling, a continuous regression, and a 5%-trimmed metric).
- `gen_dumps.sh` — regenerates the per-seed dumps from the shipped Ginnie
  design under the recommended config (needs `build-release/decision_mesh`).
- `RESULTS_*.txt` — captured console output.

## Reproduce
```bash
cmake -S core -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release -j
experiments/shrinkage_prototype/gen_dumps.sh                 # writes /tmp/reshrink_repro/seed{11,13,21}_*
python3 experiments/shrinkage_prototype/robust.py            # cross-seed table
```

## Method
Each admitted free coefficient contributes a raw surplus `delta_raw`, recovered
from the dump by undoing the stored MAP shrinkage
(`delta_pooled = delta_raw / (1 + lambda*s)`), with data variance `s = sigma_sq`.
The prior model is `delta_raw ~ N(0, tau_cell^2 + s)`. Schemes differ only in
how `tau` is indexed. Scoring splits the admitted coefficients into fit/score
halves 200x; each scheme's `tau` is fit on the fit half and its EB marginal
log-likelihood is evaluated on the held-out half. This rewards a prior that
predicts the *dispersion of unseen coefficients*, and controls for the extra
parameters the cell schemes spend.

## Findings (4 meshes: shipped baseline n=183; re-run seeds 11/13/21 n~1300)
1. **Validation.** Offline per-depth ML `tau` reproduces the C++ stored
   shrinkage per depth to ~0.01 (except tiny-n depth 10). The offline model is
   faithful.
2. **The exchangeability spread is real.** Within a single depth, K_eff spans
   ~3 orders of magnitude (e.g. depth 6: 0.007 -> 69), yet the current rule
   assigns one `tau` to all of them.
3. **Cells help — but only with enough coefficients per cell.**
   Held-out EB log-lik delta vs per-depth:
   - shipped baseline (n=183): fixed bins **-67**, pooled **-47** (cells overfit)
   - seed 11 (n=1262): **+11 / +11**
   - seed 13 (n=1466): **+13 / +18**
   - seed 21 (n=1300): **+1 / +6**
   On realistic mesh sizes the cells win by +6 to +18 nats; on the tiny mesh
   they crater. The cross-dump mean is dragged negative entirely by the one
   anomalously small shipped dump.
4. **Partial pooling and quantile bins are the robust choice.** Fixed bins (B)
   are brittle (worst on small n); tertile bins with pooling toward the depth
   mean (Q) capture most of the large-mesh gain while limiting small-mesh
   damage. The naive continuous log-tau regression (R) is misspecified and
   discarded.

## Implication for the redesign
Move `tau` to *(depth x K_eff quantile)* cells **with partial pooling toward
the depth mean**, and **fall back to per-depth when a cell is sparse** — the
small-mesh failure mode is real and must be guarded. Knobs: bin count, pooling
strength `kappa`, sparse-cell threshold.

## Caveat
The metric is held-out *coefficient* calibration (necessary, not sufficient).
The production target is held-out *surface deviance/NLL*, which depends on the
shrunk heights and needs a real run (C++ port, or a surface-reconstruction
re-score emitting `run_obs.csv`). This prototype sizes and de-risks the change;
it does not replace that measurement.
