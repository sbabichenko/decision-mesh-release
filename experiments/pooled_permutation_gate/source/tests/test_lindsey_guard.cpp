#include "lfdr.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

double normal_cdf(double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

double normal_quantile(double probability) {
    double lo = -6.0;
    double hi = 6.0;
    for (int iteration = 0; iteration < 80; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        if (normal_cdf(mid) < probability) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

bool all_equal_one(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return value == 1.0;
    });
}

}  // namespace

int main() {
    // A smooth deterministic normal sample should produce a normalized,
    // concave Lindsey density and finite local-FDR values.
    std::vector<double> normal_scores;
    normal_scores.reserve(2000);
    for (int index = 0; index < 2000; ++index) {
        const double probability = (index + 0.5) / 2000.0;
        normal_scores.push_back(normal_quantile(probability));
    }
    const LfdrResult valid = lindsey_lfdr(normal_scores, true);
    if (!valid.used_lindsey || valid.status != LindseyStatus::valid ||
        std::fabs(valid.fitted_mass - 1.0) > 0.02 ||
        !(valid.central_curvature < 0.0) ||
        !std::all_of(valid.lfdr.begin(), valid.lfdr.end(), [](double value) {
            return std::isfinite(value) && value >= 0.0 && value <= 1.0;
        })) {
        std::cerr << "valid Lindsey fit failed guard invariants\n";
        return 1;
    }

    // A degenerate histogram previously allowed the Poisson IRLS result to be
    // consumed even when it did not converge. It must now fail closed so the
    // caller can use its declared BH fallback.
    const LfdrResult invalid = lindsey_lfdr(std::vector<double>(100, 0.0), true);
    if (invalid.used_lindsey || invalid.status == LindseyStatus::valid ||
        !all_equal_one(invalid.lfdr)) {
        std::cerr << "invalid Lindsey fit did not fail closed\n";
        return 2;
    }

    return 0;
}
