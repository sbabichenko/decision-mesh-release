# Approach catalog

Every approach uses the same exact coordinate solver and fixed production safeguards unless explicitly stated: piecewise-linear solver on, hierarchical fit/map on, 12-logit height limit, four-logit IRLS proposal cap, q=.10, and deterministic execution.

## Gate approaches

### `legacy`: analytic composite baseline

Uses the legacy candidate coordinate and analytic variance correction, followed by the existing empirical-null/local-FDR machinery. This remains the current baseline because it behaved best on the matched homogeneous certified null, despite being less objective-consistent than the exact score.

### `b3`: exact signed-root gain with free Lindsey empirical null

Scores a candidate by `sign(delta_hat)*sqrt(2*penalized_objective_gain)`. At every adaptive round, Lindsey's Poisson histogram regression estimates the marginal score density and central matching estimates the empirical null. The new guard rejects nonconverged, nonfinite, badly normalized, or nonconcave fits and uses direct BH for that round.

### `b3fixed`: constrained empirical-null scale

Same exact score and Lindsey density, but constrains the empirical-null standard deviation to [1,1.4]. This tests whether free scale estimation is causing over-selection. It helped the clean null but did not help Gaussian or heavy-tail pool effects.

### `b3bh`: direct normal BH

Skips Lindsey and converts exact scores to two-sided standard-normal p-values before BH. This is a deliberately simple comparator, not a theoretically justified null for every stage. Its poor Gaussian result quantifies how costly the universal N(0,1) assumption can be.

### `histmix`: matched heavy-tail histogram null

Uses a deterministic histogram/null law matched to the generating 90/10 two-normal mixture in each size stratum. It measures the attainable gain from knowing the nuisance law rather than estimating a fresh empirical null from adaptively selected scores.

## Intervention approaches

### Adaptive stopping grid

The same frozen data and seed are stopped after `(stages, maximum rounds per stage)` equal to `(1,1)`, `(1,3)`, `(3,1)`, and `(3,3)`. These are diagnostic boundaries, not production fits. Pairing by seed isolates the incremental effect of within-stage recursion and repeated stage-wise reuse.

### Parent-uncertainty scale

The propagated parent variance term is multiplied by 0, .25, or .50. The .25 value is the reviewed default. This is a path-changing intervention: later candidate families differ after early decisions, so results measure total algorithmic effect, not merely a static denominator change.

### Exposure equalization

The Gaussian complete-null experiment is repeated after replacing each row exposure by 394, approximately the original mean, while retaining locations and the same declared random-effect law. This separates heterogeneous information weights from overdispersion.

### Shared-dependence stress

Independent Gaussian effects are replaced by a common-plus-idiosyncratic construction within randomized 40-row blocks while holding marginal variance at .35. The comparison measures sensitivity to ignored dependence at rho=.01 and .05.

### Smooth-truth decomposition

A known quadratic certified interior is embedded within a smooth global surface. Admissions can be legitimate at shallow levels, so the primary error metric is truth-projected nonquadratic texture; raw admissions and face counts are treated as workload/geometry diagnostics.

## Lindsey correctness guard

A Lindsey fit is usable only if Poisson IRLS converges, coefficients/densities are finite, fitted mass is in [0.98,1.02], and central log-density curvature is negative. Failure invokes BH and is logged. The regression test contains a valid deterministic normal sample and a degenerate histogram that must fail closed.

## What the comparisons identify

These interventions estimate total changes in the realized adaptive algorithm. They do not form an additive ANOVA decomposition: mechanisms interact, and changing an early decision changes all descendant candidates. The comparisons are therefore paired causal interventions on the algorithm, not claims that the reported effects sum to a single total.
