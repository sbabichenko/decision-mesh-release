# Predeclared null / benchmark instruments (v3, v4, and the battery)

These are the exact CSVs behind every pinned certification count in
the program. They are data, not generators: treat them as frozen
instruments. Provenance metadata where available: *_truth.npz /
*_meta.npy / null_truth_beta.npy.

## Files
- ginnie_certified_null_v3.csv: the primary certified null (real
  Ginnie covariate/exposure structure, null signal). PIN: 2
  deep-in-window under the legacy production config.
- ginnie_null_v4_pipeline.csv: pipeline-matched v4 null. PIN: 3-7
  (multi-seed range). ginnie_null_v4.csv is the pre-pipeline
  variant (byte-identical size; use _pipeline for acceptance).
- *_off.csv: offset-column variants (5th csv column), used for the
  offset-neutrality end-to-end checks. Offset is retired on the
  mainline; keep for regression only.
- cos_g{0,10,32,69}.csv + *_truth.npz: cosine signal family, g0 is
  the null member; truth files give the injected surface for power/
  recovery checks.
- confound.csv, knull_k*.csv, chirp meta: the remaining battery
  members (chirp inputs are regenerable from chirp_meta.npy).

## The deep-in-window counting protocol (the acceptance metric)
Counts CERTIFIED-FALSE deep structure: faces at depth >= 8 whose
centroid lies in the predeclared window WALA in [105, 175], WAC in
[3.5, 5.5], on a null run. Verbatim reference implementation
(coordinates rescaled by the REAL data ranges, not the null file's):

    m = pd.read_csv(dump_prefix + '_mesh.csv')
    d0 = pd.read_csv('real_pools_all.csv')
    ar = 0.5*np.abs((m.x1-m.x0)*(m.y2-m.y0)-(m.x2-m.x0)*(m.y1-m.y0))
    dep = np.round(np.log2(0.5/np.maximum(ar,1e-12))).astype(int)
    cx = (m.x0+m.x1+m.x2)/3*(d0.wala.max()-d0.wala.min())+d0.wala.min()
    cy = (m.y0+m.y1+m.y2)/3*(d0.wac.max()-d0.wac.min())+d0.wac.min()
    count = int(((dep>=8)&(cx>=105)&(cx<=175)&(cy>=3.5)&(cy<=5.5)).sum())

## How the pins were produced (context, not a bit-exact target)
Pins v3=2 / v4=3-7 were derived under the LEGACY clamped estimator
with the composite gate: DMESH_GATE_EXTRA_VAR=0.35, q=0.10, CORE
envs (SIGMA0_MIN=0.3, SIGMA0_MAX=6, PI0_MIN=0.05,
TAU_FLOOR_LOGIT_SD=0.005), DMESH_QUAD_PRIOR=1. Your B3 calibration
is a different solver and currency: the acceptance is that YOUR
re-derived thresholds put v3 at ~2-class and v4 inside 3-7 on these
exact files, i.e. the same false-certification scale on the same
instruments, not byte-equality with a retired configuration.
Run pattern reference: DMESH_DATA=<null.csv> DMESH_SPLIT=1 <binary>
0 0.10 1 24 7 (plus your stack envs).

## Reminder: v3/v4 is necessary, not sufficient
Completion of phase B requires all four acceptance criteria from
HANDOFF_SESSION2.md: v3 ~ 2 and v4 in 3-7, t_split 0/20, base <= -6.9,
NLL within 0.01 of 2.9347. The three regression splits reported so
far are not the battery.
