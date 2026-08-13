# TAU_FLOOR bypass: source read, fix, and A/B measurement (2026-08-12)

Resolves the bug lead queued in the August 6 session (writeup
sec:aug6-ladder; `experiments/aug6_session/SESSION_NOTES_2026-08-06.md`):
"`TAU_FLOOR` bypass at depths 8-9 (`lambda_v` = 1e8 wholesale vs the
floor-implied cap of 4; d10 shows the floor working). Source read queued."

## Mechanism (source read)

`lambda_v` is exactly `1 / tau_for_depth(depth)` (`core/vertex.cpp`,
`compute_eb_params`), and every `tau_for_depth` read applies
`max(tau, tau_floor_sq)`. With `DMESH_TAU_FLOOR_LOGIT_SD=0.005` the floor is
`1e4 * 0.005^2 = 0.25` centilogit^2, i.e. the cap `lambda_v <= 4`. There is
no per-depth bypass. What happened instead:

The final post-selection re-shrinkage pass (`DMESH_FINAL_RESHRINK=1`,
`core/pipeline.cpp`) re-estimated tau per depth from raw surpluses with the
truncation correction and then **replaced the floor itself**:
`mesh->tau_floor_sq = 1e-8;`. That opened two weld channels:

1. **Direct.** A depth whose truncation-corrected moment equation reads
   `chi <= target(0)` gets tau = 0, floored at 1e-8, i.e. `lambda_v = 1e8`:
   every admitted coefficient at that depth is hard-welded in one stroke,
   regardless of its individual information. This is the wholesale
   depth-8/9 weld measured in the August 6 session.
2. **Persistent.** The 1e-8 floor outlives the pass. The post-reshrink
   height sweeps call `maybe_recompute_tau_sq()` on every commit interval
   (`core/mesh.cpp`), so any later in-loop tau recompute that reads zero at
   a depth also welds it -- even when the reshrink pass itself read a
   positive tau there. Measured below (e.g. SMM seed 107 depth 6: reshrink
   read tau = 48.5, final dump shows `lambda_v = 1e8`).

The August 6 "d10 shows the floor working (lambda = 4)" observation is the
improvement-only clamp `t = min(t, tau_for_depth(d))` passing through an
in-loop value that had already been floored at 0.25 -- not the floor
holding at read time.

## Fix

`core/pipeline.cpp`: the reshrink pass keeps its honest re-estimation
(truncation correction, improvement-only clamp, `4^{-dd}` decay borrowing)
but no longer lowers `tau_floor_sq`. A tau = 0 moment read now shrinks to
the configured floor (`lambda_v = 4` at the documented sd = 0.005) instead
of welding, and the floor stays in force for all subsequent recomputes.
`DMESH_RESHRINK_NO_FLOOR=1` restores the historical behavior exactly
(verified byte-identical against a pre-fix build; with
`DMESH_FINAL_RESHRINK` unset the binaries are byte-identical too, so runs
without the reshrink pass are untouched).

## A/B measurement

Harness and raw results: `experiments/tau_floor_fix_aug12/`. Arms per
(dataset, seed): `legacy` = fix + `DMESH_RESHRINK_NO_FLOOR=1` (historical
behavior, bit-verified), `fixed` = fix, `off` = no reshrink. Binary args
`0 0.10 1 24 <seed>`.

### June 2026 SMM design (recommended one-law config + reshrink, 11 seeds: 7, 101-110)

Held-out mean deviance/pool:

| arm    | mean   | sd     |
|--------|--------|--------|
| legacy | 1.1661 | 0.0541 |
| fixed  | 1.1552 | 0.0386 |
| off    | 1.1545 | --     |

The floor bound on 8 of 11 seeds. Fixed improves 7, ties 3 (bit-identical
runs where no depth read below the floor), and is +0.001 on one (seed 108,
noise level). Largest single effect, seed 103: the depth-6 **workhorse
tier** (40 coefficients, median |surplus| 43 centilogits) read tau = 0
after the truncation correction and was welded wholesale by the legacy
pass; held-out deviance 1.292 -> 1.205. This reproduces the August 6 weld
disease on shipped data, and at a coarse signal-bearing depth, not only in
thin deep tiers. Note also that with the fix the reshrink pass is
approximately deviance-neutral against no-reshrink (1.1552 vs 1.1545);
most of the legacy pass's held-out cost was the weld, not the honest
shrinkage.

### April proxy, benchmark config + reshrink (seeds 7, 101-105)

Floor binds on 2 of 6 seeds (small deep tiers at depths 9-10); held-out
deviance changes are +-0.001. Negligible either way.

### Shipped flat true-null (`tests/data/flat_p78_clean_seed0.csv`, 5 split seeds)

Three splits admit nothing (arms identical). Two splits admit ~126
depth-4 coefficients whose tau reads 0.2 (< floor); legacy welds all of
them, fixed caps at lambda = 4. Held-out deviance: fixed equals the
no-reshrink arm and is never worse than legacy (mean delta -0.0002). The
fix does not degrade null-regime behavior: at these scales the floored
prior (sd = 0.005 logits) still shrinks essentially to the parent
interpolation, it just does so through a finite, priceable lambda instead
of an unpriceable weld.

### Generated SMM-geometry flat nulls (3 datasets x 2 split seeds)

All arm pairs bit-identical (every tau read positive). The instrument's
absolute levels are off-calibration (crude k ~ Binomial(n0, 0.852%)
generator; see `experiments/tau_floor_fix_aug12/`), so it is recorded for
completeness, not as evidence.

## Relation to the August 6 EB audit

This fix removes the accidental component of the binary trust ladder: no
tier can be welded to `lambda_v = 1e8` by a single under-supported moment
read, and every admitted coefficient keeps a finite, priceable stiffness.
It does not implement the release rule or the spike-and-slab ladder design
(writeup sec:aug6-release / sec:aug6-ladder); those remain the principled
successors. The 2020-regime June->July fits that motivated the bug lead
use monthly SFPS inputs not shipped in this archive, so the wholesale
depth-8/9 weld is reproduced here on the June 2026 SMM design instead
(8 of 11 seeds, including a full workhorse-tier weld at seed 103).
