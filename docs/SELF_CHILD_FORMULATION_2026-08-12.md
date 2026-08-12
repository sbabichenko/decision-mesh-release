# Self-child formulation (DMESH_SELF_CHILD=1) - 2026-08-12

## Motivation

In the shipped architecture a vertex height `h_v` plays three roles at once:
(1) the best coarse approximant over its admission-time footprint, (2) the
final surface's nodal value, and (3) the prior anchor for every descendant's
shrinkage null. These optima coincide only where the truth is locally linear -
exactly where refinement would not have happened. The joint Gauss-Seidel solve
lands `h_v` on a compromise, and the decision layer (gate, certificates, EB
ladder, attribution) reads that compromise as if it were each role separately.
Documented symptoms: the sawtooth lineages ("teeth, not the bite"), the
gerrymander biography, the 2.2x joint-fit-correlation se optimism, and the
silent, un-gated repricing of coarse heights ("coarse scales are not
thresholded").

The self-child formulation makes the vertex at round k+1 a child of itself at
round k. The surface becomes a telescoping sum of per-round correction
increments; a coefficient is tested exactly once, against a fixed offset, and
never silently re-fit afterward.

## Semantics

- After the stage's pre-round base fit (the level-0 model), every coefficient
  freezes (`sc_ref = height`, `mesh.sc_frozen = true`).
- Each gate round the candidate family is uniform: (a) inactive edge midpoints
  (new geometry), (b) active conformity-only welds (releases), and (c) active
  free coefficients (self-children, i.e. gated re-adjustments). Every
  candidate's null is its CURRENT rendered value ("this round changes nothing
  here"); for self-children the EB prior also centers there, so the priced and
  shrunk object is the delta. New midpoints/welds keep the configured prior
  target (e.g. the quadratic inherited-child center under DMESH_QUAD_PRIOR),
  which is a frozen constant at scoring time.
- Admissions run through the unchanged gate (B3/lindsey/BH). Commits are
  sequential with the existing live re-profile; each accepted vertex thaws
  (`sc_thawed = true`) with its delta reference (`sc_ref`) and commit-time
  prior center (`prior_mean`) recorded.
- After the round's commits, one joint polish solves ONLY the round's
  admissions (exact tempered binomial objective, diagonal prior at the
  commit-time centers - the coupled-stencil machinery does not apply because
  prior centers are frozen constants). Then everything re-freezes and each
  admitted vertex books `delta_pooled = height - sc_ref` (the realized
  increment - the per-depth tau ladder therefore learns from round increments,
  the natural exchangeable population in this formulation).
- No stage-level joint re-solves, no HIER_FIT=2 final refit, no final
  retirement, no FINAL_RESHRINK: the frozen-increment surface is the model.
- A REJECTED self-delta changes nothing: the existing coefficient keeps its
  freedom, height, and historical admission record.

## What this buys (design intent)

- Certificates never mutate: the admission-time statistic describes the
  shipped increment forever; FDR accounting is genuinely sequential.
- Coarse repricing is charged and logged: re-adjustment of an existing height
  is an explicit admission, not a free silent re-fit (the gerrymander becomes
  impossible as a silent event; its repair channel is the default).
- Welds and unproposed candidates and existing coefficients are the same
  object (a constraint priceable through the gate), completing the release
  rule of sec:aug6.
- The coupled hierarchical penalty dissolves (prior anchors are constants
  within any solve).

## Known costs and v1 caveats

- Greedy vs joint: coarse structure discovered late must be expressed as many
  small deltas; level-k estimation noise must be re-priced (or locked in) at
  later rounds. This is the accepted trade; the A/B measures it.
- Tau tiers remain keyed by geometric depth, not by round; parent-variance
  term (0.25*sp) in the candidate null retained even though parents are frozen
  (conservative).
- Root vertices (no parent_edge) freeze after the base fit and are never
  candidates (also true of the baseline solver's coordinate set).
- Dormant coefficients are not offered as self-children.
- The b1 fixed-point audit reports nonzero joint gradients by design (the
  frozen surface is not the joint MAP).
- The empirical-null population now mixes fresh candidates and self-deltas;
  central matching absorbs this by construction, but the mixture should be
  watched in the B3 calibration lines.

## Implementation map

- `core/units.h`: `self_child_mode()` env check.
- `core/vertex.h`: `sc_thawed`, `sc_ref` fields.
- `core/mesh.h`: `sc_frozen` flag.
- `core/pipeline.cpp`:
  - freeze point after pre-round sweeps;
  - `pl_solver::run`: coordinate set restricted to thawed when frozen;
    diagonal prior at `prior_mean`; coupled stencil off; starved-root
    override off;
  - legacy sweep branch: same restriction;
  - frontier: active free vertices added as delta candidates;
  - scoring: null/prior/legacy-center overrides for self-children;
  - commit: self-children allowed; rejection carve-out; thaw bookkeeping;
  - post-admission polish + re-freeze with realized-delta booking;
  - final passes (HIER_FIT=2 refit, retirement, FINAL_RESHRINK, TAX_DEST,
    final sweeps) skipped.

Baseline behavior with the flag unset is bit-identical (all changes are
guarded).
