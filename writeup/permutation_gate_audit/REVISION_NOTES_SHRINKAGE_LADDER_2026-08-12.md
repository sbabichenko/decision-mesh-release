# Revision note: sec:aug12 (shrinkage ladder, exchangeability, two-currency verdict)

Adds Section~\ref{sec:aug12} (`sec_aug12_shrinkage.tex`, \input after the
Aug 6 session) recording the August 12 empirical-Bayes implementation work
against the sec:aug6-ladder audit:

- resolves the queued TAU_FLOOR weld bug (floor overwrite in the reshrink
  pass; ledger entry 12);
- the reshrink shrinkage-ladder family (cells, keffcont, eb, ebscale,
  twogroups) and the exchangeability principle (K_eff is data support, not
  signal: it belongs in the likelihood via tau/(tau+v), not in the prior
  grouping). eb (selection-corrected MMLE + hierarchical borrow) is the
  recommended reshrink estimator for the point-prediction currency;
- honest per-pool variance (size de-bucketing of the one-law strata);
- the two-currency verdict: the ladder improves held-out deviance (best
  1.1502 vs 1.1545 naked on SMM) but is inert in the marginal NLL (full
  April Ginnie, 120,595 pools: all estimators within 2e-4), because
  re-shrinkage edits heights at frozen topology and cannot change admission;
- the in-loop prior swap (DMESH_INLOOP_EB) that does reach admission
  over-admits and regresses both currencies (SMM seed 7: mesh 348->454,
  dev 1.189->1.435, NLL 0.9074->0.9387), identifying the null fraction pi on
  the full candidate ensemble as the missing regulariser and upgrading the
  unified spike-and-slab from preference to requirement.

State-page addendum and ledger updated. Supporting harness and raw records:
`experiments/tau_floor_fix_aug12/`; standalone notes:
`docs/TAU_FLOOR_RESHRINK_FIX_2026-08-12.md`,
`docs/HIERARCHICAL_SHRINKAGE_LADDER_2026-08-12.md`.

PDF not regenerated in this environment (no LaTeX toolchain); the .tex
sources are the current record.
