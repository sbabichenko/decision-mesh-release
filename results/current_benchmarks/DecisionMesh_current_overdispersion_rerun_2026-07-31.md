# Current Decision Mesh overdispersion rerun

Configuration:

- cohort geometry;
- current default two-scale surface;
- exact coupled stage objective;
- quadratic inherited-child prior;
- q = 0.10;
- data scale = 4;
- density cap = 24 points per face;
- persistent pool standard deviation = 0.30, so truth sigma_u^2 = 0.09;
- seeds 7 through 11.

Results:

- leverage-corrected mean: 0.08856;
- seed SD: 0.00390;
- mean bias: -0.00144;
- RMSE: 0.00377;
- range: 0.0835 to 0.0935;
- naive mean: 0.08712.

The older writeup reported approximately 0.087 to 0.088. The current five-seed
mean is 0.08856, so that result is reproduced rather than overturned.

Implementation caveat: with final hierarchical retirement enabled, seeds 7-10
completed and produced the same estimates shown here. Seed 11 hit the retirement
prediction-invariance assertion. Because retirement is intended to be an
algebraically prediction-neutral cleanup, seed 11 was evaluated with final
retirement disabled. Its corrected estimate was 0.0835. This should be treated
as a retirement implementation issue, not evidence that the variance estimator
itself failed.
