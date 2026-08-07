# Hierarchical support lifecycle (Phase B1)

The old cleanup classified support from a vertex's shrinking nodal hat. That
criterion is invalid once stages are fitted in hierarchical coordinates: an
ancestor can have an empty nodal star while its coefficient still changes
training predictions through constrained descendants.

The current implementation uses one shared hierarchical-column builder
(`hierarchical_design.{h,cpp}`) for the exact stage solver and the fixed-topology
support census.

## Three separate support notions

- **Structural support:** the raw squared norm of the hierarchical column on
  training observations. This alone decides retirement.
- **Likelihood information:** the current weighted curvature of that column.
  It controls fitting precision but not whether the coordinate exists.
- **Certification support:** effective pool/issuer replication. It belongs to
  the gate, not to retirement.

## Adaptive-stage dormancy

Before each EB depth-scale update, current support diagnostics are recomputed
on hierarchical columns. A free non-root coefficient is dormant only when its
full training column is empty. The same pass updates its hierarchical
`current_fit_xTx`, point count, `sigma_sq`, and effective support. Thus a
nodally empty ancestor supported through descendants remains awake and enters
the EB census in the coordinate system actually fitted by the exact solver.

## Fixed-topology retirement

After the gate closes:

1. Build the training hierarchical design using the explicit even/odd split,
   never `mesh.weights` as a proxy for train membership.
2. Identify free non-root coefficients whose raw column norm is zero to
   numerical precision.
3. Retire only the deepest empty depth in each batch.
4. Materialize every constrained midpoint at zero surplus in birth order.
5. Assert that the batch changes training predictions by at most `1e-7`
   centi-logits before any re-optimization.
6. Rebuild the design and repeat to a fixed point.
7. Re-estimate EB scales excluding retired coefficients and run one final exact
   fixed-topology solve.

Deepest-first rebuilding matters because retiring a descendant removes a stop
in the ancestral column; an ancestor that initially appeared empty can regain
support and must not be retired.

## Candidate profiling

Exact candidate optimization no longer requires six rows. Any nonempty
training column with positive EB precision has a finite one-dimensional MAP;
replication requirements are deferred to gate certification. Remaining
unscoreable candidates are reported by cause. On the original split all 12
were true no-training-column candidates, with zero tempered-out candidates.

## Original-split regression (July 28, 2026)

- held-out mixture NLL: `2.9347`
- WALA 0-2 law-marginal dot / fitted surface: `-7.95 / -7.91`
- retired: 67 coefficients in 8 deepest-first batches
- nodally empty but hierarchically supported coefficients woken: 12
- maximum immediate training-prediction change in every batch: exactly 0
- exact candidate profiles: 18,236 scored, 12 no-training-column,
  0 tempered-out, 0 negative gains

Regression seeds remain stable:

- seed 202: NLL `2.9525`
- seed 404: NLL `2.9435`
