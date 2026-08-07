# Phase B0 post-B1 scaffolding-removal results

Date: 2026-07-28

## Tested configuration

The completed B1 line was exercised with:

```text
DMESH_PL_SOLVER=1
DMESH_HIER_FIT=2
DMESH_HIER_MAP=1
DMESH_SPLIT=1
```

The three legacy safeguards were disabled one at a time relative to the pinned
configuration (`height limit=12 logits`, `IRLS step cap=4 logits`, and the
saturation-triggered exact fallback enabled in the pre-cleanup source).

A safeguard is considered removable only if the original split and regression
seeds 202/404 preserve the fitted path to numerical tolerance and the 20-split
`t_split.py` battery remains stable.

## Original split

| configuration | law-mixture NLL | fitted WALA 0-2 | active vertices | free coefficients | gate admissions | max |height| |
|---|---:|---:|---:|---:|---:|---:|
| pinned B1 | 2.9335 | -6.96 | 2696 | 1702 | 1640 | 11.679 |
| height limit off | 2.9448 | -7.02 | 2601 | 1845 | 1788 | 16.575 |
| IRLS step cap off | 2.9409 | -6.84 | 2639 | 1782 | 1723 | 11.684 |
| saturation-exact off | 2.9335 | -6.96 | 2696 | 1702 | 1640 | 11.679 |

The law-marginal WALA 0-2 reference is -7.95.

### Height guard verdict: active, retain

The 12-logit guard rejected 106 precommit proposals on the original split.
Disabling it changed the selected topology, increased gate admissions by 148,
allowed fitted coefficients to reach 16.57 logits, and moved held-out NLL by
+0.0113. It is not inert and is retained until the gate statistic is migrated
off the legacy working scale.

### IRLS step cap verdict: active, retain

Disabling the four-logit working-step cap changed the selected topology and
admissions, moved held-out NLL by +0.0074, and weakened the original-split
boundary from -6.96 to -6.84, failing the <= -6.9 boundary acceptance pin. It
is not inert and is retained until the gate no longer consumes IRLS working
responses.

### Saturation-exact verdict: regression-inert, remove

With the saturation trigger disabled, the original split and seeds 202/404 had
identical discrete hierarchy state: vertex IDs, coordinates, parents, depths,
active/free/retired/dormant flags, and gate-admission flags all matched. Across
those three runs:

- maximum fitted-probability difference: 1.01e-10;
- maximum final height difference: 1.01e-7 centi-logits;
- active/free/admission counts: identical;
- reported predictive metrics: identical to printed precision.

The tiny numerical differences came from an intermediate scratch proposal that
is superseded by B1's mandatory live exact re-profile before commitment. The
fallback no longer changes the estimator and has been deleted, together with
its now-dead saturation-pause state and `DMESH_ADMIT_PAUSE` hook.

## Regression seeds

| split | pinned NLL | height-off NLL | step-off NLL | fitted base: pinned / height-off / step-off |
|---|---:|---:|---:|---:|
| 202 | 2.9498 | 2.9499 | 2.9496 | -6.99 / -6.83 / -7.02 |
| 404 | 2.9450 | 2.9452 | 2.9388 | -6.72 / -6.73 / -6.76 |

The seed-specific boundary values are regression diagnostics; the <= -6.9
acceptance pin applies to the original split.

## Twenty-split battery after removing saturation-exact

With the saturation-triggered fallback disabled:

```text
t_split: eruptions 0/20 | mean 2.9449 | sd 0.0054 | worst 2.9537
PASS
```

The observed standard deviation is below the pinned 0.0067 reference and well
inside the predeclared two-times sanity ceiling of 0.0134. No split exceeded the
3.2 eruption threshold.

## B0 decision

- Delete saturation-triggered exact fallback: **yes**.
- Delete associated saturation admission-pause state: **yes**.
- Delete 12-logit height guard: **no; still causally active in the legacy gate path**.
- Delete four-logit IRLS step cap: **no; still causally active in the legacy gate path**.

B0 therefore removes one dead mechanism. The two remaining safeguards are not
part of the exact committed-coordinate solver; they remain because candidate
ranking and local-FDR calibration still use the legacy IRLS/Wald gate. Their
removal should be retested after B2/B3 migrate that gate.
