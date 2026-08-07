# Coupled hierarchical-penalty correction - July 30, 2026

## Problem

The stage optimizer stores canonical vertex heights but penalizes hierarchical
surpluses. A parent-height move changes the surplus residuals of free descendants
whose inherited parents contain that height. The earlier exact-stage implementation
included only the visited coordinate's own Gaussian term, so its coordinate fixed
point need not be stationary for the declared penalized objective.

For free-surplus vector `delta = A h`,

    P(h) = 0.5 * (A h - m)' Lambda (A h - m),
    grad P(h) = A' Lambda (A h - m).

## Implementation

`core/pipeline.cpp` now:

1. builds the sparse derivative of every penalized free surplus with respect to
   every canonical free height;
2. uses all affected penalty terms in each one-dimensional safeguarded Newton
   profile;
3. incrementally updates the affected prior residuals after a committed
   Gauss-Seidel move; and
4. optionally audits the full objective gradient and boundary KKT conditions when
   `DMESH_COUPLED_FIXED_POINT_AUDIT=1`.

The stage-initial exact fit is also enabled by removing the early return at refresh
epoch zero. The change is active for the default linear parent-midpoint prior. The
existing quadratic-prior and ghost-boundary modes retain their previous specialized
treatment.

## Validation

The clean release build passed the Lindsey fail-closed regression test.
The validated Ginnie run bundled under `docs/coupled_penalty_validation/` reports:

| Variant | Marginal NLL |
|---|---:|
| Release baseline | 2.9410 |
| Stage-initial exact fit, diagonal-only prior | 2.9520 |
| Stage-initial exact fit, coupled prior | 2.9443 |

On the final fit, 1,256 nonroot coordinates were audited. The 1,255 interior
coordinates had maximum absolute gradient `3.65e-6` nats per logit and RMS
`1.64e-7`; the single boundary coordinate satisfied KKT.

This is retained as a correctness repair. It improves materially over the
nonstationary exact-stage variant but does not beat the release baseline overall.
