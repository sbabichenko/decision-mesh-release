# Real Ginnie Mae fit after removing WALA 0

## What was done

Two checks were performed across ten matched deterministic held-out splits:

1. The existing full-data fits were rescored after excluding WALA 0.
2. WALA 0 was removed from the input data and both gates were refit from
   scratch.

The filtered input has 120,135 pools, down from
120,595; 460 WALA-0 pools were
removed.

## Existing fits, WALA 0 removed only from scoring

| Gate | Held-out deviance/pool | X2/information |
|---|---:|---:|
| Shipped gate | 6.5527 | 0.15078 |
| Permutation BH | 6.8673 | 0.15650 |

The deviance gap falls to
0.3145 per pool.

## Full refit after removing WALA 0

| Gate | Held-out deviance/pool | X2/information | Admissions | Faces | Max depth |
|---|---:|---:|---:|---:|---:|
| Shipped gate | 6.4968 | 0.15237 | 118.0 | 729.8 | 8.7 |
| Permutation BH | 6.8108 | 0.15857 | 136.1 | 885.9 | 6.1 |

The permutation gate remains worse in all ten splits. Its mean paired
deviance increase is 0.3140, with a
split-bootstrap interval of [0.2328,
0.3970].

Removing WALA 0 reduces the original mean deviance gap from about 0.5615 to
0.3141, so the suspect WALA-0 observations explain roughly 44% of the apparent
loss, but not all of it.

The shipped gate's mean admitted coefficient count falls from 187.7 on the
full data to 118.0 after removing WALA 0. The
permutation count changes much less, from 143.1 to
136.1. This supports the view that the shipped gate
spent substantial model capacity fitting the WALA-0 boundary feature.

The largest remaining gap is at WALA 1--5, especially WALA 1 and 2, with
additional smaller losses at high WALA and high coupons.
