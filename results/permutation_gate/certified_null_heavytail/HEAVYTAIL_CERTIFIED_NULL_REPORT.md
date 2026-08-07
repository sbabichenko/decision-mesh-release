# Heavy-tailed corrected certified-null experiment

## Design

This experiment keeps the corrected certified-null spatial truth and Ginnie
WALA/WAC exposure geometry. The truth is exactly quadratic inside WALA
[105,175] and WAC [3.5,5.5], while genuine smooth signal remains outside that
window.

The Gaussian pool effect was replaced by a centered 90/10 two-normal mixture:

- total variance: 0.35;
- core variance: 0.184211;
- wide variance: 1.842105;
- wide/core variance ratio: 10;
- theoretical kurtosis: 9.058.

Twenty deterministic seeds were generated. These are new mechanism-equivalent
stress tests, not members of the earlier frozen acceptance set.

## Results across 20 seeds

| Configuration | Total admissions | Window admissions | Seeds with window admissions | Deep window admissions | Deep faces | Window nonquadratic RMS | Outside centered RMSE | Outside correlation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Shipped gate, q=.10 | 165.60 | 13.40 | 18/20 | 2.70 | 115.85 | 11.13 | 16.41 | 0.9870 |
| Permutation BH, K_eff≥.5, q=.10 | 4.70 | 0.00 | 0/20 | 0.00 | 0.00 | 3.38 | 14.73 | 0.9905 |
| Permutation BH, K_eff≥4, q=.10 | 3.80 | 0.00 | 0/20 | 0.00 | 0.00 | 3.39 | 14.23 | 0.9913 |
| Permutation BH, K_eff≥4, q=.15 | 7.75 | 0.00 | 0/20 | 0.00 | 0.00 | 3.37 | 14.12 | 0.9914 |
| Permutation BH, K_eff≥4, q=.20 | 16.45 | 0.25 | 4/20 | 0.00 | 3.30 | 3.91 | 14.00 | 0.9914 |

## Main findings

The heavier-tailed certified null exposes a real failure of the shipped gate.
It admitted structure inside the exactly quadratic window in 18 of 20 seeds,
made deep admissions in 7 seeds, and produced an average of 115.85 deep faces
inside that window. Its average nonquadratic RMS was 11.13 centi-log-odds, with
a maximum of 42.62.

Permutation BH with K_eff >= 4 and q=.10 or q=.15 made no admissions and no
deep faces inside the certified window in any of the 20 seeds. Both also
improved the fit to the genuine signal outside the window relative to the
shipped gate.

The q=.15 setting is the strongest clean compromise in this sweep:

- zero certified-window admissions in 20/20 seeds;
- zero deep certified-window faces in 20/20 seeds;
- 7.75 total admissions on average;
- certified-window nonquadratic RMS 3.37;
- outside centered RMSE 14.12, versus 16.41 for the shipped gate.

Increasing q to .20 lowers outside RMSE slightly further, to 14.00, but it
reintroduces certified-window admissions in 4 of 20 seeds and raises the
window texture metric to 3.91. The paired difference in outside RMSE between
q=.20 and q=.15 is small and uncertain, while the texture penalty is positive
on average. Therefore q=.15 is the preferable operating point on this battery.

## Caveats

The permutation implementation still uses 32 deterministic pooled permutation
draws per candidate rather than the proposed double-saddlepoint tail
calculation. The heavy-tailed datasets were generated specifically for this
test, so they are reproducible but were not part of the original frozen
acceptance package.
