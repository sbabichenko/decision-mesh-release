# Revision notes - pooled-permutation gate and Ginnie outcome audit

Date: August 1, 2026

Base source: `DecisionMesh_writeup_posterior_pool_effects_2026-08-01.tex`.

## Main changes

1. Corrected the April Ginnie response description.
   - The single snapshot does not observe SMM.
   - `round(Original Aggregate / AOLS) - current loan count` is now described as a constructed loan-count attrition proxy, not exact cumulative terminations.
   - Removed the old real-data figure whose embedded labels called the response cumulative termination probability.
   - Recast the older WALA-boundary section as historical solver forensics on an unavailable snapshot.

2. Added the WALA-zero sensitivity audit.
   - Removed 460 WALA-zero pools and refit both gates on ten matched locally balanced splits.
   - Added the held-out comparison table and paired-deviance figure.
   - Recorded that WALA zero explains about 44% of the original permutation-gate performance gap, while the remaining gap is concentrated first at WALA 1-5.

3. Added the experimental pooled-permutation gate.
   - Documented exact first-three permutation moments, 32 deterministic pooled draws per candidate, direct BH, and the commit-time effective-support check.
   - Clearly labels the implementation as experimental and distinct from the proposed deterministic double-saddlepoint tail calculation.

4. Added the new batteries.
   - 34-run frozen complete-null battery.
   - Five-seed Gaussian corrected certified null.
   - Twenty-seed heavy-tailed certified null.
   - Ten-split real Ginnie comparison before and after the WALA-zero exclusion.

5. Updated status, discussion, limitations, and next steps.
   - The permutation gate now has strong tested null control but insufficient real-data power.
   - It is not adopted as the release default.
   - The next gate experiments are depth-aware support and a coarse-to-fine hybrid, followed by enumeration-validated double-saddlepoint tails.

6. Removed four placeholder figure panels whose source graphics were never supplied. Their numerical conclusions remain in the prose.

## Deliverables

- `DecisionMesh_writeup_permutation_gate_audit_2026-08-01.tex`
- `DecisionMesh_writeup_permutation_gate_audit_2026-08-01.pdf`
- Updated supporting figures and summary CSVs.
