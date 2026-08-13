# tau-floor reshrink fix A/B (2026-08-12)

Findings and tables: `docs/TAU_FLOOR_RESHRINK_FIX_2026-08-12.md`.

- `sweep.py`: runs the three arms (legacy = `DMESH_RESHRINK_NO_FLOOR=1`,
  fixed, off = no reshrink) over the June 2026 SMM design (seeds 7,
  101-110), the April proxy (seeds 7, 101-105), the shipped flat true-null
  (5 split seeds), and optional generated SMM-geometry nulls. Repo root is
  hardcoded at the top. Dumps are summarized (held-out deviance, admitted
  count, per-depth lambda_v ladder, [reshrink] tau lines) and deleted.
- `analyze.py`: renders `results_all.json` into the per-run and
  divergent-run tables (`analysis.txt`).
- Generated nulls: k ~ Binomial(n0, 0.00852) on the real SMM design
  geometry, `random.Random(1000 + i)`, inverse-CDF sampler; regenerate with
  the snippet in the git history or any binomial sampler. Off-calibration
  as an absolute instrument; arms were bit-identical on it.

Equivalence checks (both byte-identical on all dump files): modified
binary vs pre-fix build with `DMESH_FINAL_RESHRINK` unset (frozen lindsey
benchmark, seed 7), and legacy arm vs pre-fix build with reshrink on (SMM
design, seed 7).

## Ladder variants (2026-08-12, follow-on)

Findings: `docs/HIERARCHICAL_SHRINKAGE_LADDER_2026-08-12.md`.

- `sweep2.py`: 6-arm ladder sweep (off, depth, cells, cells_rawkeff,
  twogroups, twogroups_free) over SMM 11 seeds + April + shipped null.
  Uses a snapshotted binary (`dm_snapshot2`) so a concurrent rebuild
  cannot swap the executable mid-run.
- `sweep_kc.py`: appends the threshold-free `keffcont` arm to
  `results_ladder.json` (adding a mode leaves the other arms numerically
  identical, so they are not re-run).
- `analyze2.py`: renders the per-seed table and two-groups EM diagnostics.
- `results_ladder.json` / `analysis_ladder.txt`: raw records and rendered
  tables for all seven arms.

## Principled EB modes (2026-08-12, follow-on)

- `sweep_eb.py`: appends the `eb`, `ebcell`, and `eb_free` arms to
  `results_ladder.json` (adding modes leaves existing arms identical). Uses
  a snapshotted binary. `analysis_ladder.txt` / `analyze2.py` render all
  arms; the recommended arm is `eb` (selection-corrected MMLE tau with a
  hierarchical global-borrow, per-vertex law-aware posterior).

## Scale prior + honest v_i (2026-08-12, follow-on)

- `sweep_scale2.py`: crosses {eb, ebscale} x {stratum v_i, smooth POINT_VAR}
  over the SMM design, plus the aggressive (DMESH_EB_FREE) ebscale arms.
- `gen_pointvar.py`: writes `pointvar_smooth_smm.txt`, a smooth log-log
  interpolation of the 4-stratum pool variance in pool size (de-bucketing),
  for DMESH_POINT_VAR.
- `results_scale2.json` / `analysis_scale.txt`: records and rendered table.
Findings and the recommended default (eb) plus the planned EB tau-multiplier
next step are in `docs/HIERARCHICAL_SHRINKAGE_LADDER_2026-08-12.md`.
