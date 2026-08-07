# Real Ginnie Mae test: pooled-permutation BH

## Configuration

Experimental gate:
- pooled-permutation direct BH;
- q = 0.15;
- minimum effective support K_eff = 4;
- 32 deterministic permutation draws per candidate;
- three stages and up to 80 gate rounds per stage.

Comparator:
- current shipped Lindsey gate at q = 0.10.

The comparison uses ten matched deterministic train/held-out splits of the same
120,595-pool Ginnie Mae cross-section. These are split-stability replicates,
not ten independent datasets.

## Ten-split results

| Metric | Shipped gate | Permutation BH | Change |
|---|---:|---:|---:|
| Held-out deviance per pool | 6.5830 | 7.1445 | +8.5% |
| Held-out X2 per information | 0.15205 | 0.16943 | +11.4% |
| Active admitted coefficients | 187.7 | 143.1 | -23.8% |
| Final faces | 1002.2 | 927.2 | -7.5% |
| Maximum active depth | 10.1 | 6.8 | -32.7% |
| Runtime | 26.0 s | 14.6 s | -43.7% |

Permutation BH had higher held-out deviance in all 10 matched splits. The mean
paired increase was 0.5615, with a split-bootstrap interval
of [0.4450, 0.6797]. It had lower X2 in one split
and higher X2 in nine; the mean increase was 0.01738.

## Where the loss occurs

About 89% of the aggregate held-out deviance loss is concentrated at WALA below
30 months. A similar share is concentrated in WAC 5--7 percent. The gate is
therefore not losing broad low-order structure uniformly; it is failing mainly
in the young-pool, mid-to-high-coupon region.

The most visible structural change is depth: mean maximum active depth falls
from 10.1 to 6.8. This supports the
interpretation that the K_eff floor plus permutation BH suppresses localized
deep refinements that are useful on the real cross-section.

## Assessment

The q=.15, K_eff>=4 rule remains excellent on the heavy-tailed certified null,
but it should not replace the shipped real-data gate in its present form. It is
too conservative for localized real Ginnie structure.

The next modification should preserve candidate-specific permutation
calibration while relaxing the hard depth/support bottleneck. Plausible tests
are a depth-aware K_eff threshold, a weaker commit-time support rule, or a
hybrid in which permutation BH governs early/coarse admissions and calibrated
lfdr governs replicated descendants.
