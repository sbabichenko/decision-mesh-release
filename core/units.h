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
