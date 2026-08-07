#pragma once

#include <cstdint>
#include <vector>

struct Vertex;

// Experimental pooled-permutation candidate calibration.
//
// Each star is standardized by its own exact finite-population permutation
// moments. Null draws are then pooled across stars of similar row count. This
// removes the per-star 1/choose(m,r) resolution floor while keeping the null
// scale candidate-specific.
struct PermutationStar {
    bool ok = false;
    int rows = 0;
    int unique_pools = 0;
    int size_band = 0;
    double raw_z = 0.0;
    double observed_z = 0.0;
    double skewness = 0.0;
    std::vector<double> null_z;
};

struct PooledPermutationResult {
    std::vector<double> p_values;
    std::vector<double> gaussian_scores;
    std::vector<int> reference_sizes;
};

[[nodiscard]] int permutation_size_band(int rows);

[[nodiscard]] PermutationStar build_permutation_star(
    Vertex& candidate,
    double null_height_units,
    int null_draws,
    std::uint64_t seed,
    bool allow_repeated_pool_rows = false);

[[nodiscard]] PooledPermutationResult pool_permutation_null(
    const std::vector<PermutationStar>& stars,
    int minimum_reference_size = 512);

[[nodiscard]] double inverse_standard_normal(double probability);
