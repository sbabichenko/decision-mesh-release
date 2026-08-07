# Decision Mesh current model audit — 2026-07-31

## Comparable 12-seed scorecard

All progression numbers below use the same seeds:
101, 102, 103, 104, 105, 106, 107, 109, 110, 111, 112, 114.

- Linear local-split mean overall NLL: 2.720933
- Quadratic child mean overall NLL: 2.715170
- Preferred nested crossed-cover mean overall NLL: 2.710123
- Quarter-precision control mean overall NLL: 2.709934

The preferred nested crossed-cover model improves overall NLL by
0.010809 relative to the matched linear baseline.
The quarter-precision control improves it by 0.010999,
but is not adopted because its precision multiplier is not yet replaced by a
joint conditional-information EB calculation.

## Current status

- Shipped July-30 code: exact coupled objective and locally balanced split;
  the geometry-generated correction is not included.
- Adopted development fix: quadratic inherited child target.
- Positive frozen-topology development: graph-patch variance sharing and nested
  crossed geometric covers.
- Experimental only: quarter interaction precision.
- Still open: exact-score false-discovery calibration, heavy-tail adaptive-null
  false branches, full 20-seed validation of the quadratic/geometric branch,
  and joint mesh-cover Schur-complement EB.
