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
