#pragma once

#include <cstdlib>

// Mesh heights are stored in centi-log-odds for numerical compatibility with
// the original implementation. Keep conversions explicit at module boundaries.
inline constexpr double kHeightScale = 100.0;
inline constexpr double kVarianceScale = kHeightScale * kHeightScale;

inline constexpr double to_log_odds(double scaled_height) noexcept {
    return scaled_height / kHeightScale;
}
inline constexpr double from_log_odds(double log_odds) noexcept {
    return log_odds * kHeightScale;
}

// DMESH_SELF_CHILD=1: sequential frozen-increment ("self-child") mode.  After
// the stage's pre-round base fit every coefficient freezes; each gate round
// prices per-vertex DELTAS (new midpoints, weld releases, and re-adjustments
// of existing coefficients alike) against the frozen surface through the same
// gate, jointly fits only the round's admissions, and re-freezes.  See
// docs/SELF_CHILD_FORMULATION_2026-08-12.md.
inline bool self_child_mode() {
    static const bool enabled = std::getenv("DMESH_SELF_CHILD") != nullptr;
    return enabled;
}

// DMESH_SC_JOINT: decouple the self-child DECISION architecture (gated
// per-vertex delta candidates) from the frozen-increment VALUE semantics.
// With SC_JOINT set, the enlarged candidate family and delta gating remain,
// but heights are re-solved jointly every round exactly as in the baseline
// (no freezing, final passes re-enabled).  1 = self-delta priors centered at
// "no change"; 2 = hierarchical prior centers (degeneracy check: predicted
// to reduce to the baseline, since a coordinate at its joint optimum has
// zero profile gain against the same objective).
inline int self_child_joint_level() {
    static const int level = [] {
        const char* e = std::getenv("DMESH_SC_JOINT");
        return e ? std::atoi(e) : 0;
    }();
    return level;
}

// True when self-child mode runs with frozen-increment value semantics
// (the original formulation).  All freezing behavior keys off this.
inline bool self_child_frozen() {
    return self_child_mode() && self_child_joint_level() == 0;
}
