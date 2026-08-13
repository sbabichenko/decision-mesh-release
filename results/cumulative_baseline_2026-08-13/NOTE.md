# Cumulative-implied baseline experiment - 2026-08-13 (MEASURED, negative)

Question: can the April cumulative book, fitted as a surface F(WALA,WAC),
serve as the baseline (offset) for the June monthly SMM fit, with the mesh
fitting only the correction?

Construction: fit F on data/ginnie_design.csv (seed-7 lindsey benchmark
surface dump); isotonic-regress F along WALA per WAC line (the raw fitted F
is non-monotone at sub-grid scale - a 1-month difference puts the MEDIAN
implied hazard at the clamp floor); derive the implied monthly hazard from
the 6-month survival ratio h = 1-(S(w)/S(w-6))^{1/6}; recenter n-weighted
to the June pooled logit (removing the proxy's ~1.6-logit level error, so
only SHAPE is imported); feed as DMESH_OFFSET5 to the June fit.

Result (seed 7, one-law config, held-out):

| arm | dev | NLL |
|---|---|---|
| no offset (baseline solver) | 1.189 | 0.9242 |
| half-strength offset | 1.231 | - |
| full offset | 1.453 | 1.0203 |
| full offset + certified self-child | 1.429 | 1.0156 |

Monotone in dose: the cumulative-implied shape is noise on net at every
strength. The mesh built a mean-|g| 1.80-logit correction (vs 0.51 without
the offset) and still could not neutralize it.

Three stacked mechanisms, each sufficient:
1. Response mismatch: the April cumulative is the retracted loan-count
   attrition proxy (pooled 0.783 - not terminations); its shape is
   dominated by whatever the proxy actually counts.
2. Age-period-cohort confounding: age-differencing a single cross-sectional
   cumulative surface yields an age-cohort mixture, not the current-period
   hazard; cohort composition drift is the same order as the age effect.
3. Mechanical: cumulative surfaces are differentiable only where monotone;
   16% of pools (5,511) sat on floored non-monotone cells, and those errors
   are pool-scattered - invisible to any smooth correction surface.

Reading: this measures, from a new direction, why the writeup separates
eras ("numbers from three eras are not cross-comparable") and why the
two-month panel is the standing instrument: WITHIN-cohort differencing
(two cumulative snapshots of the same pools) is the valid hazard estimator;
cross-sectional age-differencing is not. If two cumulative tapes a month
apart are available, the per-pool difference IS the monthly count
attrition, and the whole certified pipeline applies directly with no
offset construction at all.

Builder script: build_offset.py in this directory.
