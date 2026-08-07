# Crossed geometric-cover empirical-Bayes experiment

## Question

Can the logic of overlapping group hierarchies be translated to Decision Mesh at the function level, rather than by attaching additional parents to split-tree surpluses?

## Baseline

The baseline is the frozen-topology quadratic-child Decision Mesh with geometry-patch EB variance scaling. The adaptive topology and all original mesh coefficients remain unchanged.

## Construction

For each fitted triangulation:

1. Build the weighted active-mesh graph from triangle edges.
2. Generate two graph-geometric covers by farthest-point center insertion, using different initial centers.
3. Give each observation soft membership in its three nearest centers under graph distance.
4. Center each cover into identifiable contrast coordinates.
5. Define the interaction basis as the tensor product of the two centered cover memberships.

The residual correction is

```text
g(x) = X_A(x) a + X_B(x) b + X_AB(x) c.
```

The three coefficient blocks have independent Gaussian EB scales. Variance components are selected on training pools only using a Laplace/IRLS approximation to the same heavy-tailed marginal pool likelihood used for evaluation. The interaction scale may collapse essentially to zero.

The selected specification uses 12 centers in each cover and three nonzero memberships per observation. It has about 131-143 effective coefficients depending on rank.

## Seed-7 screen

With 12 centers, relative to the geometry-patch baseline NLL 2.711093:

| Model | Overall NLL | Middle-aged NLL | Seasoned NLL |
|---|---:|---:|---:|
| Parent covers only | 2.707293 | 3.044284 | 2.368035 |
| Interaction only | 2.706460 | 3.043216 | 2.367628 |
| Parents + interaction | **2.706454** | **3.043113** | **2.367628** |

Changing the second cover's initial center among `(0.5,0.5)`, `(1,0)`, `(0,1)`, and `(1,1)` changed overall NLL only from approximately 2.70645 to 2.70649.

## Matched 12-seed results

Changes are relative to the quadratic + geometry-patch baseline. Negative is better.

### Parent covers only

| Region | Mean NLL change | Improved seeds |
|---|---:|---:|
| Overall | -0.003813 | 12 / 12 |
| Ramp | -0.001317 | 11 / 12 |
| Middle-aged | -0.001436 | 12 / 12 |
| Seasoned | -0.008634 | 12 / 12 |

### Full crossed model

| Region | Mean NLL change | Improved seeds |
|---|---:|---:|
| Overall | **-0.004100** | **12 / 12** |
| Ramp | -0.001866 | 11 / 12 |
| Middle-aged | **-0.001907** | **12 / 12** |
| Seasoned | **-0.008703** | **12 / 12** |

### Increment from interactions beyond parent covers

| Region | Mean NLL change | Improved seeds |
|---|---:|---:|
| Overall | -0.000287 | 11 / 12 |
| Ramp | -0.000549 | 10 / 12 |
| Middle-aged | -0.000472 | 11 / 12 |
| Seasoned | -0.000069 | 8 / 12 |

The interaction EB scale exceeded 0.001 in 10 of 12 splits; it collapsed to the numerical floor in two splits. Thus the interaction is useful but not universally required.

## Stability

Correction fields fitted on the 12 different training splits have mean pairwise correlation 0.902, with range 0.844 to 0.963. The RMS of the cross-seed mean correction is about 3.34 times the RMS cross-seed fluctuation, indicating a stable omitted function-level component rather than split-specific noise.

## Interpretation

This is the first hierarchy experiment with a practically meaningful predictive gain. Most of the improvement comes from the two broad overlapping geometric covers. Their intersections add a smaller but repeatable improvement, especially in the middle-aged band.

The result supports the blog-post logic in a more principled form:

- parent geometries are directly estimated function-level effects;
- intersections are explicit crossed interactions;
- overlap is handled jointly in the likelihood rather than through an ad hoc effective-sample correction;
- split-tree surpluses are not treated as comparable group estimates.

## Caveats and next step

This is currently an external frozen-topology correction. It has not yet been integrated into candidate admission or the mesh's joint C++ objective. Before production integration, the correction should be represented as a sparse partition-of-unity basis in the final refit, with the three EB variance components estimated jointly. Adaptation should remain unchanged until a separate calibration study shows that using the cover model inside the gate is stable.
