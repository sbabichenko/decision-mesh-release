# Coupled hierarchical-penalty fix

The coordinate solver previously priced only the visited coefficient's own
Gaussian prior. A parent move also changes the surpluses—and therefore the
penalties—of free descendants. The patch constructs the sparse derivative of
every penalized surplus with respect to every canonical free coordinate and
updates those residuals during Gauss-Seidel sweeps.

## Fixed-point validation

On the final Ginnie fit:

- 1,256 nonroot coordinates were audited.
- 1,255 were interior coordinates.
- Maximum interior full-objective gradient:
  3.65e-6 nats per logit.
- Interior gradient RMS:
  1.64e-7 nats per logit.
- One coordinate was at the -12-logit height cap.
- Boundary KKT violations: zero.

Thus the patched solver is stationary for the coupled penalized objective,
up to numerical tolerance and the declared height constraint.

## Ginnie result

The approved marginal NLL changes as follows:

- Release baseline: 2.9410
- Stage-initial exact fix alone: 2.9520
- Stage fix plus coupled penalties: 2.9443

The coupled correction recovers about 70% of the NLL deterioration introduced
by enabling the skipped stage fit:

    (2.9520 - 2.9443) / (2.9520 - 2.9410) = 0.70

It does not beat the release baseline overall. It improves ramp NLL and
ordinary held-out deviance, while the seasoned-region NLL is worse.

## Synthetic controls

Against the stage-initial exact fix, the coupled penalty changes the three-seed
mean centered RMSE only in the fifth decimal:

- clean flat null: 0.004730 -> 0.004730
- Gaussian-pool null: 0.124234 -> 0.124230
- smooth truth: 0.214302 -> 0.214254
- corrected certified truth: 0.220047 -> 0.220008

The topology and admission counts are unchanged in these shallow three-round
controls. The omission matters primarily in the much deeper Ginnie hierarchy.

## Status

This is a mathematically correct prototype for the default linear
parent-midpoint prior. It intentionally falls back to the old treatment when
quadratic-prior or ghost-boundary modes are enabled; those priors need their
own derivative construction. It should not yet replace the release baseline
without retuning the adaptive gate against marginal NLL.
