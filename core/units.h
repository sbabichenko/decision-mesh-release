#pragma once

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
