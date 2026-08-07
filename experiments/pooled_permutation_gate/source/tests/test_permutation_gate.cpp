#include "permutation_gate.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace {

bool nearly_equal(double a, double b, double tolerance) {
    return std::fabs(a - b) <= tolerance;
}

}  // namespace

int main() {
    // Inverse-normal round trips at common quantiles.
    if (!nearly_equal(inverse_standard_normal(0.5), 0.0, 1e-10) ||
        !nearly_equal(inverse_standard_normal(0.975), 1.95996398454, 2e-7) ||
        !nearly_equal(inverse_standard_normal(0.001), -3.09023230617, 2e-7)) {
        std::fprintf(stderr, "inverse normal regression failed\n");
        return 1;
    }

    // A pooled reference should remove the per-star resolution floor. Each
    // candidate contributes only 16 draws, but 400 candidates give a reference
    // of 6400 draws and therefore p-values well below 1/16.
    std::mt19937_64 rng(12345);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::vector<PermutationStar> stars(400);
    for (auto& star : stars) {
        star.ok = true;
        star.rows = 24;
        star.size_band = permutation_size_band(star.rows);
        star.observed_z = normal(rng);
        for (int draw = 0; draw < 16; ++draw) star.null_z.push_back(normal(rng));
    }
    // Plant one strong candidate.
    stars[0].observed_z = 4.0;
    const auto pooled = pool_permutation_null(stars, 512);
    if (pooled.p_values.size() != stars.size() ||
        pooled.reference_sizes[0] < 6000 ||
        !(pooled.p_values[0] < 0.01) ||
        !(pooled.gaussian_scores[0] > 2.5)) {
        std::fprintf(stderr, "pooled calibration regression failed: p %.8g ref %d z %.8g\n",
                     pooled.p_values[0], pooled.reference_sizes[0],
                     pooled.gaussian_scores[0]);
        return 1;
    }

    // Under a matched Gaussian null, pooled p-values should be approximately
    // uniform. This is intentionally loose: it is a deterministic smoke test,
    // not a stochastic certification.
    int below_05 = 0;
    for (double p : pooled.p_values) below_05 += p <= 0.05;
    if (below_05 < 8 || below_05 > 35) {
        std::fprintf(stderr, "null-size smoke test failed: %d/400\n", below_05);
        return 1;
    }

    return 0;
}
