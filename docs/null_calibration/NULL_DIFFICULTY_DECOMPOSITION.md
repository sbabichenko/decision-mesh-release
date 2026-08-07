# Quantitative null-difficulty decomposition

Date: 2026-07-28

## Executive result

The Lindsey failure has been fixed as a fail-closed software invariant. The patched code accepts a Lindsey local-FDR fit only when Poisson IRLS converges, all coefficients and density values are finite, the fitted density integrates to 0.98--1.02, and the central log-density is concave. Otherwise it uses direct BH and records the failure reason.

The controlled battery shows that the remaining problem is not one undifferentiated "hard null." Each number below is tied to an explicit frozen scenario and an explicit algorithmic intervention shipped with this package. The mechanisms have very different measured sizes:

1. **Catastrophic numerical Lindsey failure:** large but rare. On the five-seed corrected homogeneous certified-null battery, the guard reduced mean total admissions from 644.2 to 233.6 (64% reduction), window admissions from 53.4 to 14.2, deep admitted coefficients from 7.0 to zero, and nonquadratic texture from 13.62 to 6.34 centi-logits.
2. **Repeated adaptive testing:** the largest universal residual effect. On a flat p=0.78 binomial null, one stage/one round admitted 0.0 coefficients across three seeds; three stages/three rounds admitted 9.67 on average and admitted something in every seed.
3. **Pool-effect law:** material but secondary to adaptation in early rounds. At 3x3, the free-Lindsey B3 gate averaged 9.00 admissions on the heavy-tail null. A matched histogram law reduced this to 1.67; the legacy analytic gate produced 0.67.
4. **Parent-estimation variance:** highly consequential when other overdispersion is absent. On the clean p=0.78 null, removing the parent term increased admissions from 9.67 to 65.67. Under Gaussian/heavy pool effects, the incremental effect was small and nonmonotone because pool-noise variance already dominated the candidate standard error.
5. **Exposure heterogeneity:** affected fitted texture more than admission count. Equalizing every row to 394 loans preserved total exposure and changed 3x3 admissions only from 9.80 to 8.67, while reducing mean certified-window nonquadratic texture from 2.49 to 1.01 centi-logits.
6. **Smooth low-order truth and geometry:** generated many legitimate shallow admissions. The smooth quadratic-window Gaussian benchmark rose from 19.3 admissions after one pass to 112.7 after 3x3, while window nonquadratic texture fell from 5.45 to 3.29. This confirms that raw face/admission counts on a smooth non-flat benchmark cannot be interpreted as candidate-level false discoveries.


## Experimental design and traceability

The package now contains the exact input scenarios and all code paths used to obtain the tables. `SCENARIO_CATALOG.md` defines every data-generating intervention; `APPROACH_CATALOG.md` defines every gate and ablation; `manifests/frozen_run_manifest.csv` maps each completed run to its input file and configuration. The exact controlled-null CSVs are shipped under `experiments/controlled_null_decomposition/frozen_inputs/`, and the exact five corrected-HCN inputs are shipped under `experiments/corrected_homogeneous_certified_null/frozen_inputs/`.

The comparisons are **paired algorithmic interventions**, not an additive variance decomposition. Holding a seed and frozen input fixed, one mechanism or stopping rule is changed and the entire adaptive path is rerun. Because early admissions alter future candidate sets, effect sizes include downstream amplification and should not be added together.

### Scenario families

- Complete flat nulls at p=.50, .78, and .97 isolate baseline probability and saturation.
- Gaussian and two-normal-mixture row effects, each with total logit variance .35, isolate overdispersion and tail law.
- Randomized 40-row common-shock blocks at rho=.01 and .05 stress ignored dependence while preserving marginal variance.
- Equalized exposure n=394 isolates heterogeneous information weights.
- A known smooth surface with an exactly quadratic certified interior separates legitimate shallow refinement from invented high-frequency texture.
- The corrected homogeneous certified null supplies the before/after Lindsey bug test on a harder, non-flat global surface.

### Approach families

- Legacy analytic composite gate.
- Exact signed-root gain with free Lindsey empirical-null calibration.
- The same score with constrained empirical-null scale.
- Direct standard-normal BH.
- A deterministic heavy-tail histogram null matched to the generating law.
- Parent-variance scales 0, .25, and .50.
- Predeclared adaptive stopping points (1x1, 1x3, 3x1, and 3x3).

See the catalogs for exact formulas, environment variables, seeds, and interpretation.

## 1. Lindsey bug: exact mechanism and measured effect

The old code could consume a nonconverged Poisson histogram fit. In the known seed-4 failure, the alleged density integrated to roughly 4,323 rather than 1. This made nearly all local FDR values approximately zero and admitted all 988 candidates in one round.

The guard is now explicit and tested. The current code reproduces the repaired full seed-4 result: 174 total gate admissions, 2 window admissions, 0 deep window admissions, and 41 deep window faces. A CTest regression checks both a valid deterministic normal sample and a degenerate histogram that must fail closed.

Across the controlled 3x3 B3 runs, guarded fallback rounds admitted **zero candidates in every condition**. All ordinary residual admissions occurred on numerically valid Lindsey rounds. Thus the guard removes catastrophic tail risk but does not by itself calibrate the adaptive procedure.

## 2. Adaptive amplification

Three paired seeds were evaluated at four predeclared stopping points:

| Null | 1 stage x 1 round | 1 stage x 3 rounds | 3 stages x 1 round | 3 stages x 3 rounds |
|---|---:|---:|---:|---:|
| Flat clean p=.78 | 0.00 | 0.33 | 4.00 | 9.67 |
| Flat Gaussian pool effect | 1.00 | 2.00 | 6.00 | 8.33 |
| Flat heavy-tail pool effect | 1.00 | 2.33 | 3.67 | 9.00 |

The stage repetition is at least as important as recursion within a stage. On the clean null, increasing one stage from one to three rounds added only 0.33 admissions on average, while repeating a single round across three stages added 4.00. The full 3x3 loop added 9.67 relative to one pass.

This points directly toward separating candidate construction from verification, rather than attempting only a more elaborate within-round null density.

## 3. Gate-law and empirical-null comparisons

Three-seed, 3x3 diagnostic-stop means:

| Null | Legacy analytic | B3 free Lindsey | B3 fixed scale box | Direct normal BH | Matched histogram |
|---|---:|---:|---:|---:|---:|
| Flat clean p=.78 | 0.33 | 9.67 | 3.67 | 6.00 | -- |
| Flat Gaussian | 2.00 | 8.33 | 8.33 | 28.33 | -- |
| Flat heavy-tail | 0.67 | 9.00 | 9.00 | 13.00 | 1.67 |

Free empirical-null scaling is responsible for part of the clean-null over-selection: fixing sigma0 to [1,1.4] reduced clean admissions from 9.67 to 3.67. It did not improve Gaussian or heavy-tail runs. Direct normal BH was particularly poor under Gaussian pool effects, so replacing Lindsey with naive BH is not a viable general solution.

## 4. Other nuisance interventions

### Baseline probability / saturation

Five-seed B3 3x3 means were 7.8 admissions at p=.50, 7.6 at p=.78, and 4.0 at p=.97. Saturation is therefore not the dominant null-admission mechanism in this guarded exact-score configuration, although it still affects coefficient magnitudes and numerical conditioning.

### Artificial shared dependence

With balanced 40-row blocks and the same marginal Gaussian variance, the three-seed 3x3 mean was 8.33 admissions under independence, 7.00 at rho=.01, and 11.67 at rho=.05. The direction at rho=.05 is adverse but the small sample is noisy. This is a stress result, not an issuer estimate; the frozen files contain no issuer identifiers.

### Parent uncertainty

| Null | Parent scale 0 | Parent scale .25 | Parent scale .50 |
|---|---:|---:|---:|
| Clean p=.78 | 65.67 | 9.67 | 5.67 |
| Gaussian | 8.67 | 8.33 | 11.33 |
| Heavy-tail | 10.33 | 9.00 | 9.00 |

The current .25 correction should be retained. The clean-null ablation shows that deleting it can produce an eruption even when the candidate's own binomial variance is well behaved. The .50 result is not a sufficient basis for retuning: under pool effects it was neutral or adverse due to nonlinear path changes.

## 5. Decision implications

1. **Merge the Lindsey guard and regression test.** This is a correctness repair, not a tuning choice.
2. **Do not treat guarded B3 as calibrated.** Numerically valid Lindsey rounds still create ordinary false admissions.
3. **Prioritize adaptive separation.** The measured stage/retest effect is larger and more universal than heavy tails, exposure heterogeneity, or moderate artificial correlation. A fixed candidate hierarchy with grouped cross-fitting and one final structured gate is now supported by evidence rather than aesthetics.
4. **Retain parent uncertainty.** It is load-bearing on the clean null.
5. **Use a law-aware null for pool effects.** Matched heavy-tail calibration reduced admissions substantially; naive normal BH made Gaussian behavior worse.
6. **Separate statistical selections from geometry.** Track admitted hierarchical surpluses and truth-projected texture; use deep-face count only as a closure-footprint diagnostic.
7. **Acquire issuer membership before claiming dependence calibration.** The rho stress test gives a sensitivity curve, not an empirical issuer correction.

## Limits

The mechanism battery uses three to five deterministic seeds per intervention and predeclared early stopping for computational comparability. The dependence experiment uses artificial balanced blocks. These results rank mechanisms and reject several directions; they are not a formal FDR proof or a precise estimate of production error rates. Because the algorithm is adaptive, the effects are generally nonadditive and can be seed-sensitive; seed-level rows are retained in `all_run_metrics.csv`.

## Reproducibility files

The package includes both the exact frozen inputs and the code used to define, run, and analyze each intervention. See `REPRODUCIBILITY.md` for commands.


- `all_run_metrics.csv`: one row per completed run.
- `adaptivity_summary.csv` and `adaptivity_paired_effects.csv`.
- `method_summary.csv`.
- `parent_uncertainty_summary.csv`.
- `nuisance_summary.csv`.
- `lindsey_round_ledger.csv` and `lindsey_round_summary.csv`.
- `lindsey_bug_effect.csv`.
