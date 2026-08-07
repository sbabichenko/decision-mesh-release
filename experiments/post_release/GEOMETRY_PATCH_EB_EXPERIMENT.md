# Geometry-patch empirical-Bayes hierarchy experiment

## Question

Can the Decision Mesh keep its quadratic split-tree coefficient coordinate while replacing depth-only EB sharing with a hierarchy generated from the final mesh geometry?

## Model

The coefficient residual remains

```text
e_v = h_v - mu_quad(v),
```

so the fitted surface, candidate geometry, quadratic child target, and coupled penalty stencil are unchanged.

At frozen topology, the active triangulation is treated as a weighted graph. A nested binary patch tree is generated recursively:

1. Within a patch, find two far-apart vertices using shortest-path graph distance.
2. Split the patch by the nearer of the two centers.
3. Recurse until patches contain approximately the target number of free coefficients.

For each non-root patch, estimate a local prior-variance multiplier `g_P` from the admitted quadratic residuals in that patch:

```text
e_v ~ N(0, sigma_v^2 + g_P tau_depth(v)).
```

The child patch's log multiplier is penalized toward its parent patch's log multiplier. The root multiplier is fixed at one, and the final leaf multipliers are centered to have global geometric mean one. Thus the extension redistributes shrinkage over graph-local regions without changing the global slab level.

The tested fine specification used:

- target leaf size: 12 free coefficients;
- minimum admitted draws per patch: 5;
- parent log-scale penalty: 6;
- variance multiplier range: 0.35 to 2.5;
- three alternating patch-scale / exact-refit blocks;
- final-only activation, so adaptation and topology are identical to the quadratic baseline.

## Seed 7

| Region | Quadratic baseline | Geometry-patch EB | Change |
|---|---:|---:|---:|
| Overall | 2.711401 | 2.711093 | -0.000308 |
| Low-WALA ramp | 2.662843 | 2.662605 | -0.000238 |
| Middle-aged | 3.045425 | 3.045364 | -0.000061 |
| Seasoned | 2.375339 | 2.375147 | -0.000192 |

The seed-7 topology is unchanged: both runs contain the same 979 triangles with identical triangle coordinates. The coupled fixed-point audit reported no boundary KKT violations; final max/RMS gradient was about `0.00314 / 0.000244`.

## Matched 12-seed comparison

Negative changes improve held-out marginal NLL.

| Region | Mean change | SD | Improved seeds |
|---|---:|---:|---:|
| Overall | -0.000129 | 0.000092 | 11 / 12 |
| Low-WALA ramp | -0.000038 | 0.000559 | 5 / 12 |
| Middle-aged | -0.000004 | 0.000010 | 7 / 12 |
| Seasoned | -0.000043 | 0.000042 | 12 / 12 |

The only overall worsening was approximately `+0.000055`. Learned leaf scales hit the lower bound `0.35` in every split, while maximum scales ranged from roughly `1.33` to `1.78`. The fitted geometries contained 14 to 23 leaf patches.

## Interpretation

This is the first tested interacting/geometric hierarchy that gives a repeatable held-out improvement. It supports the user's concern that split-tree genealogy is not the best structure for EB borrowing: graph-local regions benefit from sharing the amount of shrinkage even when their vertices are distant in the split tree.

The gain is small and is not primarily a solution to the middle-aged WAC ridges. Its strongest consistency is in the seasoned region. It is therefore best viewed as a general regularization improvement rather than a ridge-specific basis correction.

The lower scale bound being active suggests that some graph patches prefer substantially stronger shrinkage than the depth-only hierarchy provides. Before relaxing the bound further, a larger seed battery and a calibration check should be run.

## Recommendation

Keep the quadratic child target. Retain geometry-patch EB as a promising final-refit option. It is preferable to the many-parent geometry hierarchy because it preserves the sparse quadratic coordinate system, adds little computational cost, and improves overall NLL in 11 of 12 matched splits.
