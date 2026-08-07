# Phase B2 exact signed-root gain results

Date: 2026-07-28

## Implemented statistic

For every scoreable candidate hierarchical coefficient, B2 now records

\[
r_v = \operatorname{sign}(\widehat\delta_v)\sqrt{2\Delta Q_v},\qquad
\Delta Q_v=Q_v(\widehat\delta_v)-Q_v(0),
\]

where `Q` is the exact tempered binomial objective plus the fixed Gaussian EB
penalty used by the B1 live commit.  The default gate is unchanged; the score
runs in shadow mode.  `DMESH_B2_USE_SCORE=1` is an explicitly untuned diagnostic
that feeds `r_v` into the existing lfdr implementation.

The score dump is selected with `DMESH_B2_SCORE_DUMP=<csv>`.  The old
`DMESH_B1_GAIN_DUMP` name remains as a compatibility fallback.

The dump includes the exact score, legacy score, exact displacement, effective
pool count, probability range, objective information, and a decomposition

\[
\Delta Q_v = \Delta \ell_v + \Delta \log \pi_v.
\]

## Numerical invariants

On the pinned original split (24,404 scoreable candidate profiles):

- maximum error in `signed_root_gain^2 = 2*gain_nats`: `5.82e-11`;
- maximum error in `data_gain + penalty_gain = total_gain`: `8.82e-11` nats;
- negative exact gains: `0`;
- signed-direction mismatches: `0`;
- a clean repeated run produced a byte-identical score CSV.

The shadow path leaves the B1 estimator unchanged:

- held-out law-mixture NLL: `2.9335`;
- WALA 0--2 fitted base: `-6.96` against the law-marginal `-7.95`.

Regression seeds remain unchanged in shadow mode:

| split | held-out NLL | fitted WALA 0--2 |
|---|---:|---:|
| 202 | 2.9498 | -6.99 |
| 404 | 2.9450 | -6.72 |

## Comparison with the legacy score

The exact score is not a cosmetic rescaling of the legacy Wald score.

On the original split:

- signed correlation: `0.253`;
- Spearman correlation of absolute ranking: `0.437`;
- sign disagreement: `31.0%` of candidate profiles;
- median within-round overlap of the top 10%: `30.8%`.

On seeds 202 and 404, sign disagreement rises to about `39%` and absolute-rank
Spearman correlation is about `0.21` and `0.20`, respectively.  This is the
expected nonlinear-regime effect: B1 commits were already exact, but the legacy
gate was still ranking candidates in a different currency.

### Named regime-straddle candidate

On the original split, the nearest candidate to raw `(WALA,WAC)=(0,4.52)` is
normalized `(0.015625,0.390625)`:

- legacy score: `-3.059`;
- exact signed-root gain: `-3.108`;
- exact penalized gain: `4.8292` nats;
- data gain: `5.1526` nats;
- EB penalty change: `-0.3234` nats;
- null/optimum: `-3.806 -> -5.575` logits;
- null probability range over the column: `0.00036--0.07987`;
- effective pool count: `4.34`.

The more revealing regression is seed 202: the nearby candidate has legacy
score `-24.49` but exact signed-root gain `+1.20`; the exact optimum moves upward
rather than downward.  The old score and live objective can therefore disagree
not only in magnitude but in direction.

## Complete-null diagnostic

A deterministic complete-null synthetic run (`DMESH_AMP=0`) shows that the raw
signed-root gain must not be treated as a universal N(0,1) pivot.

| stage / round regime | exact-score behavior |
|---|---|
| stage 0 | sd about `4.3--6.7`, large tails |
| stages 1--2 | center near zero, sd about `0.61--0.89` |

The stage-0 inflation has a concrete cause: stage 0 intentionally uses the
untempered binomial weights before the pool-variance estimate is available, so
pool heterogeneity enters the exact objective as apparent spatial evidence.
Later stages use the estimated pool tempering and are much closer to a
standardized null, though still not exactly unit scale.

Directly feeding the score into the old lfdr gate on this complete null caused a
large transient stage-0 mesh (`1,559` vertices), followed by final later-stage
collapse to the coarse mesh.  It did not create a final false admission, but the
transient proves that B3 needs stage/fit-state-aware calibration rather than a
global N(0,1) assumption.

## Untuned production-score experiment

Using `DMESH_B2_USE_SCORE=1` with the existing empirical-null machinery on the
original split gave:

- held-out NLL `2.9169`;
- WALA 0--2 fit `-7.00`;
- active vertices `3,379` versus `2,696` on the pinned legacy gate;
- gate-admitted final vertices `2,268` versus `1,640`.

This is predictively encouraging, but it is not a certification result.  On
seed 202 the untuned exact-score path did not finish within the five-minute
closure budget after creating a much larger adaptive workload.  The score is
therefore retained in shadow mode by default.

## Missing predeclared instruments

The supplied B1/B0 tree does not contain either certified-null data file:

- `ginnie_certified_null_v3.csv`;
- `ginnie_null_v4_pipeline.csv` (older notes call it `ginnie_null_v4.csv`).

Consequently the required v3/v4 shadow distributions and certified-window
counts could not be rerun from this bundle.  The implementation and analysis
harness are ready for those files, but B2 cannot be declared fully closed on
its predeclared instruments until they are restored.

## B2 decision

The signed-root exact penalized gain is accepted as the candidate-ranking
coordinate to carry into B3.  It is exact, deterministic, direction preserving,
and objective consistent.  It is explicitly **not** accepted as a standard
normal pivot.

B3 must address at least:

1. stage-0/pool-variance calibration;
2. fit-state and round dependence;
3. effective-replication and heavy-tail effects;
4. a pipeline-matched empirical or permutation null;
5. the selection-threshold metadata used by the EB truncation correction.

The default estimator remains the B1 legacy-calibrated gate until that work is
complete.
