# Locally balanced Ginnie holdouts and 20-seed audit

## Problem

The historical `DMESH_SPLIT=1` convention used row parity: even CSV rows trained the mesh and odd rows were held out. `ginnie_design.csv` is ordered in WALA/WAC, so neighboring high-exposure observations could be systematically separated. A large held-out miss could then appear in a region where the training half supplied no corresponding refinement signal. The split also existed only implicitly; several evaluation utilities reconstructed it from row number.

The archived 2.9443 result was produced on an older Ginnie snapshot whose run log reports pooled rate 0.7851. The CSV shipped in the release has pooled rate 0.783241. The older input is unavailable, so 2.9443 is retained only as a historical coupled-penalty validation result and is not presented as reproducible from the shipped archive.

## Implementation

For `DMESH_SPLIT_MODE=local` (now the real-data default):

1. Keep every pool intact.
2. Recursively median-partition pool centroids in WALA/WAC until a local block contains at most 24 pools.
3. Within each block, sort by exposure and form an equal-cardinality greedy partition that minimizes exposure imbalance.
4. Randomly flip the two block labels using the requested seed. Outcomes are not used.
5. Store the resulting training mask explicitly in `Dataset`, copy it into `DecisionMesh`, and use it everywhere the fitter, exact hierarchical design, fixed-point audit, output writer, and held-out evaluator need split membership.
6. Write `is_train` in `run_obs.csv`; evaluation and plotting utilities consume that field rather than reconstructing parity.

The historical behavior remains available with `DMESH_SPLIT_MODE=parity`.

## Same-file seed-7 comparison

Both runs below use the shipped 120,595-row `data/ginnie_design.csv` and the same coupled-penalty executable.

| split | marginal NLL | ramp NLL | mid NLL | seasoned NLL | faces | active vertices |
|---|---:|---:|---:|---:|---:|---:|
| row parity | 2.7413 | 2.920 | 3.045 | 2.382 | 307 | 181 |
| locally balanced | 2.7129 | 2.665 | 3.046 | 2.376 | 967 | 526 |

This is not a claim that changing the test set improves the estimator by 0.0284 NLL. It shows that the ordered parity split was a distorted evaluation/design path: under the locally representative split, the training gate sees the low-WALA structure that also appears in holdout.

## Twenty seeded local holdouts

Seeds 101-120 were run from scratch. These are repeated seeded holdouts, not an iid nonparametric bootstrap; their spread measures split sensitivity.

- Marginal NLL: mean **2.7191**, SD **0.0093**, range **2.7105-2.7543**.
- Descriptive 2.5/50/97.5 percentiles: **2.7116 / 2.7172 / 2.7441**.
- Ramp NLL: mean **2.7105**, SD **0.0836**.
- Face count: mean **976.5**, range **234-1507**.
- Low-WALA, WAC 5.5-6.25 weighted mean Pearson z: mean **+0.159**, range **-1.639 to +2.525**. The original one-sided patch does not recur systematically.

Eighteen runs form a tight main cluster: mean NLL **2.7164**, SD **0.0023**, range **2.7105-2.7193**. Seeds 101 and 107 are genuine topology-collapse outliers: they stop at 234 and 315 faces, with NLL 2.7543 and 2.7328. They are included in the primary 20-seed summary. This is a separate adaptive-gate instability exposed by the better split, not a reason to revert to row parity.

## Files

- `local_split_bootstrap/bootstrap_20seed_metrics.csv`
- `local_split_bootstrap/bootstrap_20seed_summary.csv`
- `local_split_bootstrap/split_comparison_seed7.csv`
- `local_split_bootstrap/bootstrap20_nll_by_seed.png`
- `local_split_bootstrap/ginnie_local_split_seed7_2x2.png`
- seed-101 and seed-107 logs for the two collapse cases

Reproduce with:

```bash
python3 core/bootstrap_local_split.py \
  ./build-release/decision_mesh \
  data/ginnie_design.csv \
  core/auditing/law_shared.json \
  bootstrap-output 4
```
