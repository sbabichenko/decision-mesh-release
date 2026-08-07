# Twenty-seed locally balanced Ginnie holdouts

Seeds: 101--120. Whole pools are kept intact. Each split uses recursive spatial blocks of at most 24 pools and equal-cardinality, exposure-balanced assignment within each block. A seeded block-wide label flip randomizes which side is training. Outcomes are not used to construct the split.

- Marginal NLL: mean **2.7191**, SD **0.0093**, range **[2.7105, 2.7543]**.
- Descriptive 2.5/50/97.5 percentiles: **2.7116 / 2.7172 / 2.7441**.
- Ramp NLL: mean **2.7105**, SD **0.0836**.
- Faces: mean **976.5**, range **[234, 1507]**.
- Active vertices: mean **529.0**, range **[141, 797]**.
- Internal executable runtime: mean **16.8 s**, range **[5.9, 30.5] s**.
- Low-WALA/WAC 5.5--6.25 weighted mean Pearson z: mean **+0.159**, range **[-1.639, +2.525]**.
- Topology-collapse runs with fewer than 500 faces: **2/20**. They remain in the primary summary.

The noncollapsed main cluster (faces at least 500) has mean NLL **2.7164**, SD **0.0023**, and range **[2.7105, 2.7193]**. The two collapsed runs are seeds 101 and 107. These are repeated seeded holdouts, not an iid nonparametric bootstrap; the spread measures split and adaptive-topology sensitivity.
