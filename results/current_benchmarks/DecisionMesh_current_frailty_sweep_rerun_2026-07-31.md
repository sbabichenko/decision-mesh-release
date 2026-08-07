# Current Decision Mesh frailty-variance sweep

## Configuration

- current exact coupled stage solver;
- exact candidate score and commit path;
- quadratic inherited-child prior;
- cohort geometry;
- q = 0.10;
- data scale = 4;
- density cap = 24 points per face;
- five seeds, 7 through 11, at each frailty level;
- no train/holdout split, matching the historical same-data leverage-calibration
  experiment;
- final retirement disabled because of the known quadratic-prior retirement
  invariant issue.

True sigma_u values were 0.0, 0.1, ..., 0.5, corresponding to variances
0.00, 0.01, ..., 0.25.

All 30 fits completed with finite estimates.

## Results

|   sigma_u |   truth |   corrected_mean |   corrected_sd |   corrected_bias |   corrected_relative_bias_percent |   naive_mean |   naive_relative_bias_percent |   admitted_mean |
|----------:|--------:|-----------------:|---------------:|-----------------:|----------------------------------:|-------------:|------------------------------:|----------------:|
|   0.00000 | 0.00000 |          0.00010 |        0.00000 |          0.00010 |                         nan       |      0.00010 |                     nan       |       509.60000 |
|   0.10000 | 0.01000 |          0.01008 |        0.00052 |          0.00008 |                           0.80000 |      0.00934 |                      -6.60000 |       382.00000 |
|   0.20000 | 0.04000 |          0.03940 |        0.00115 |         -0.00060 |                          -1.50000 |      0.03774 |                      -5.65000 |       256.00000 |
|   0.30000 | 0.09000 |          0.08656 |        0.00363 |         -0.00344 |                          -3.82222 |      0.08414 |                      -6.51111 |       164.40000 |
|   0.40000 | 0.16000 |          0.15166 |        0.00536 |         -0.00834 |                          -5.21250 |      0.14846 |                      -7.21250 |       123.40000 |
|   0.50000 | 0.25000 |          0.22950 |        0.01143 |         -0.02050 |                          -8.20000 |      0.22590 |                      -9.64000 |        95.80000 |

The origin-constrained calibration slopes were:

- naive: 0.9132;
- leverage-corrected: 0.9306.

## Interpretation

The leverage correction improves calibration at every positive frailty level.
It is nearly unbiased for sigma_u <= 0.2 and materially reduces downward bias
at larger frailty. It does not, however, eliminate the bias uniformly:

- at truth 0.09, corrected mean is 0.08656 (-3.8%);
- at truth 0.16, corrected mean is 0.15166 (-5.2%);
- at truth 0.25, corrected mean is 0.22950 (-8.2%).

The historical statement that the correction removes an approximately 9%
downward bias at every scale is therefore too strong for the current estimator.
A more accurate statement is that it substantially reduces same-data leverage
bias, with residual underestimation increasing at high frailty.

The number of admitted vertices falls sharply with frailty, from a mean of
509.6 at the complete null to 95.8 at sigma_u=0.5. This suggests part of the
remaining high-frailty bias may arise from interaction between robust
downweighting, reduced refinement, and variance attribution, rather than from
the trace correction alone.
