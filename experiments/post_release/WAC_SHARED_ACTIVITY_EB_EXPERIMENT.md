# WAC shared-activity EB experiment

## Model

The quadratic geometric child target is unchanged. For each admitted surplus

\[
e_v = h_v - \mu_{\mathrm{quad},v},
\]

the depth-only slab variance is multiplied by a partially pooled WAC-strip scale:

\[
e_v \sim N\!\left(0,\; \sigma_v^2 + g_{c(v)}\tau_{d(v)}^2\right).
\]

The strip scales are estimated at frozen topology from admitted training-side coefficients, with a Gaussian penalty on `log(g_c)`, normalized to have weighted geometric mean one, and capped. This redistributes shrinkage across WAC strips without changing the global slab level.

Tested specification:

- quadratic child target;
- 0.25-percentage-point WAC strips;
- minimum 4 admitted coefficients per strip;
- log-scale penalty 8;
- scale range 0.35 to 2.5;
- six alternating scale/refit blocks;
- final-only activation, so the topology is identical to the quadratic baseline.

## Results

Across the 12 completed quadratic-baseline split seeds:

| Region | Mean NLL change | SD | Seeds improved |
|---|---:|---:|---:|
| Overall | +0.000002 | 0.000066 | 5/12 |
| Low-WALA ramp | -0.000048 | 0.000649 | 6/12 |
| Middle-aged ridge band | -0.000006 | 0.000010 | 8/12 |
| Seasoned | -0.000008 | 0.000017 | 9/12 |

Negative is better. The learned scales were nontrivial: final minimum scales ranged from 0.427 to 0.609 and maximum scales from 1.341 to 1.903. Therefore the null result is not because the EB fit collapsed all scales to one.

## Interpretation

The model detects heterogeneous WAC-strip activity, but changing only the amount of shrinkage at an already selected, frozen topology does not materially improve held-out fit. The middle-aged gain is approximately six millionths of NLL per pool, far below the roughly 0.005 improvement obtained by a direct training-only WAC-by-WALA correction.

This rejects the simple hypothesis that the ridge miss is mainly caused by a depth-only slab variance that is too small on active coupon lines. More likely explanations are:

1. the relevant structure is WAC-by-WALA localized rather than a strip-wide variance effect;
2. the gate fails to admit the collective ridge basis in the first place, so final-only rescaling cannot recover it;
3. activity should be shared at the group-score or function-component level rather than through post-selection coefficient variances.

The patch is an experimental ablation and is not recommended for the production baseline.
