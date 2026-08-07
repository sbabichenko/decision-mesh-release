# Pooled-permutation gate on the corrected certified null

Five frozen seeds were run for each configuration. The truth is exactly quadratic inside WALA [105,175] and WAC [3.5,5.5], with genuine low-order signal outside that window and independent Gaussian pool effects.

All configurations produced zero gate admissions at depth 8 or greater inside the certified window.

## Mean results

| Configuration | Total admissions | Window admissions | Deep window faces | Window nonquadratic RMS | Outside centered RMSE | Outside correlation |
|---|---:|---:|---:|---:|---:|---:|
| Shipped gate, q=.10 | 110.0 | 8.2 | 73.8 | 4.69 | 13.47 | 0.9913 |
| Permutation BH, K_eff≥.5, q=.10 | 25.2 | 0.4 | 1.6 | 4.38 | 13.60 | 0.9916 |
| Permutation BH, K_eff≥4, q=.10 | 18.8 | 0.2 | 1.6 | 3.69 | 13.16 | 0.9922 |
| Permutation BH, K_eff≥.5, q=.20 | 44.4 | 0.6 | 3.2 | 4.43 | 13.70 | 0.9913 |
| Permutation BH, K_eff≥4, q=.20 | 30.0 | 0.2 | 0.8 | 4.18 | 13.18 | 0.9922 |

## Interpretation

- The current shipped gate already avoids direct deep admissions in the certified window, so this is not the catastrophic failure seen in the complete-null battery.
- Permutation BH nevertheless reduces total admissions and almost removes unnecessary certified-window refinement.
- The K_eff >= 4, q=.10 version lowers nonquadratic texture inside the quadratic window from 4.69 to 3.69 centi-log-odds.
- Outside the certified window, where real spatial signal is present, that same version slightly improves centered RMSE (13.47 to 13.16) while maintaining weighted correlation above .992.
- Raising q to .20 increases admissions without improving the five-seed fit metrics, so q=.10 is preferable on this battery.

The present implementation still uses 32 deterministic pooled permutation draws rather than the proposed double-saddlepoint tail calculation.