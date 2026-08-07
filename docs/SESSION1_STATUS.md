# Session 1 status: state semantics and containment

Completed:

1. Separate geometric presence from statistical freedom.
2. Create conformity vertices at exact zero hierarchical surplus.
3. Keep geometry-only vertices eligible for future gate promotion.
4. Validate coordinate proposals before propagating their state.
5. Mark exactly zero-support admitted coefficients dormant during adaptive topology changes.
6. Wake dormant coefficients automatically if support returns.
7. Retire still-zero-support non-root coefficients only after topology is fixed.
8. Hold near-unidentified nodal coordinates when current training rows are below 10 or effective support is below 0.5.
9. Remove the default hard height ceiling; retain `DMESH_HEIGHT_LIMIT_LOGIT` only as a diagnostic option.
10. Verify matched 12-logit, 15-logit, and unlimited real-data runs.

Final matched-run result:

| Height setting | Held-out deviance/pool | Active vertices | Active height range |
|---|---:|---:|---:|
| 12 | 12.621 | 2,514 | -8.534 to +7.134 logits |
| 15 | 12.621 | 2,514 | -8.534 to +7.134 logits |
| unlimited | 12.621 | 2,514 | -8.534 to +7.134 logits |

All three retire 40 zero-support coefficients. The 12- and 15-logit surface dumps are identical. The unlimited surface differs from the capped surfaces by at most 0.00160 logits on the evaluation grid and 0.00309 logits at observed points.

Next sessions:

- Session 2: residual-error and exact-live-versus-quadratic diagnostics.
- Session 3: true hierarchical-surplus design columns.
- Session 4: monotone exact penalized-likelihood coordinate solver.
- Session 5: likelihood-native refinement and gate recalibration.
