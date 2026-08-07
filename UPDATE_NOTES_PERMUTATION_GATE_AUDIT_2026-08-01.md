# August 1 update - pooled-permutation gate audit

## Stable production path

The production estimator in `core/` remains the adopted exact coupled-penalty,
locally balanced, quadratic inherited-child configuration. It was rebuilt from
scratch and its Lindsey regression test passed.

## Experimental branch

`experiments/pooled_permutation_gate/` contains the independently implemented
pooled-permutation branch, a patch against `core/`, unit tests, and clean build
logs. Both the Lindsey and permutation tests passed. The branch is not enabled
by default.

## Main evidence added

- 34 frozen complete-null datasets: permutation BH with
  `K_eff >= 4` made zero admissions in the three-round and full-adaptation
  batteries.
- Five-seed Gaussian corrected certified null: strong reduction in unnecessary
  certified-window refinement.
- Twenty-seed heavy-tailed certified null: `q=.15`, `K_eff >= 4` made zero
  certified-window admissions and improved outside-window fit.
- Ten matched real-Ginnie splits: the same rule underfit the real cross-section.
- WALA-zero audit: after removing 460 WALA-zero pools and refitting, mean
  held-out deviance was 6.4968 for the shipped gate and 6.8108 for permutation
  BH. WALA zero explained about 44% of the original gap, but the permutation
  gate remained worse in all ten splits.

## Outcome correction

The bundled April cross-section does not contain SMM. The development response
uses `round(Original Aggregate / AOLS) - current loan count`; this is now
reported as a constructed loan-count attrition proxy, not exact cumulative
terminations. The latest writeup removes the misleading real-data figure and
qualifies all April interpretations accordingly.

## Latest writeup

See `writeup/permutation_gate_audit/` for the updated 77-page PDF, LaTeX source,
figures, revision notes, and summary tables.

# August 5 update - exact-solver default and one-law audit

See results/aug5_pl_onelaw_audit/SESSION_NOTES_2026-08-05.md and the updated
writeup (writeup/permutation_gate_audit/DecisionMesh_writeup_aug5_audit_2026-08-05.pdf,
81pp). Headline: DMESH_PL_SOLVER is now the recommended default; the one-law
gate with the converged SMM law and the empirical-null sd floor is the
recommended admission configuration; the oscillation-era heavy-tail law is
retracted; DMESH_LAW_WEIGHTS is removed from recommended configs.
