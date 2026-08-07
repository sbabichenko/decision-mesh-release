#include "permutation_gate.h"

#include "face.h"
#include "mesh.h"
#include "observation.h"
#include "units.h"
#include "vertex.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_set>
#include <vector>

namespace {

constexpr double kTiny = 1e-14;

struct StarRows {
    std::vector<double> a;
    std::vector<double> g;
    int unique_pools = 0;
    bool repeated_pool_rows = false;
};

StarRows collect_star_rows(Vertex& candidate, double null_height_units) {
    StarRows out;
    auto faces = candidate.active ? candidate.get_faces(false)
                                  : candidate.get_sim_faces(false);
    if (faces.faces.empty() || candidate.mesh == nullptr ||
        candidate.mesh->obs_store == nullptr) {
        return out;
    }

    DecisionMesh& mesh = *candidate.mesh;
    const auto& observations = *mesh.obs_store;
    std::unordered_set<int> seen_rows;
    std::unordered_set<int> seen_pools;

    for (Face* face : faces.faces) {
        if (face == nullptr || face->coords_indices.empty()) continue;
        int slot = -1;
        for (int j = 0; j < 3; ++j) {
            if (face->vertices[j] == &candidate) {
                slot = j;
                break;
            }
        }
        if (slot < 0) continue;

        Vertex* v0 = face->vertices[0];
        Vertex* v1 = face->vertices[1];
        Vertex* v2 = face->vertices[2];
        const double determinant =
            (v1->x - v0->x) * (v2->y - v0->y) -
            (v2->x - v0->x) * (v1->y - v0->y);
        if (std::fabs(determinant) < kTiny) continue;

        double heights[3] = {v0->height, v1->height, v2->height};
        heights[slot] = null_height_units;

        for (int row : face->coords_indices) {
            if (row < 0 || row >= mesh.n_points ||
                row >= static_cast<int>(observations.size())) {
                continue;
            }
            if (!seen_rows.insert(row).second) continue;
            const double weight = mesh.weights[row];
            if (!(weight > 1e-9) || !std::isfinite(weight)) continue;

            const double x = mesh.get_x(row, 0);
            const double y = mesh.get_x(row, 1);
            const double l0 =
                ((v1->x - x) * (v2->y - y) -
                 (v2->x - x) * (v1->y - y)) /
                determinant;
            const double l1 =
                ((v2->x - x) * (v0->y - y) -
                 (v0->x - x) * (v2->y - y)) /
                determinant;
            const double l2 = 1.0 - l0 - l1;
            const double basis = slot == 0 ? l0 : slot == 1 ? l1 : l2;
            if (!(basis > 1e-8) || !std::isfinite(basis)) continue;

            const double prediction =
                l0 * heights[0] + l1 * heights[1] + l2 * heights[2];
            const double residual = mesh.values[row] - prediction;
            // Preserve the slot-specific heteroskedastic variance and permute
            // only standardized residual shape.  The original score is
            //   sum a_i w_i e_i
            // = sum [a_i w_i sqrt(v_i)] [e_i/sqrt(v_i)].
            // This is essential on the Ginnie geometry, where exposure (and
            // therefore sampling variance) is strongly spatially structured.
            double pool_variance_logit = 0.0;
            if (!mesh.gate_s2.empty()) {
                pool_variance_logit = std::max(mesh.gate_s2[row], 0.0);
            } else if (const char* extra = std::getenv("DMESH_GATE_EXTRA_VAR")) {
                const double parsed = std::atof(extra);
                if (std::isfinite(parsed) && parsed > 0.0)
                    pool_variance_logit = parsed;
            }
            const double variance_units = 1.0 / weight
                + pool_variance_logit * kVarianceScale;
            if (!(variance_units > 0.0) || !std::isfinite(variance_units)) continue;
            const double standard_deviation = std::sqrt(variance_units);
            const double slot_coefficient = basis * weight * standard_deviation;
            const double standardized_residual = residual / standard_deviation;
            if (!std::isfinite(slot_coefficient) ||
                !std::isfinite(standardized_residual)) continue;

            out.a.push_back(slot_coefficient);
            out.g.push_back(standardized_residual);
            const int pool = observations[row].pool_id;
            if (!seen_pools.insert(pool).second) out.repeated_pool_rows = true;
        }
    }
    out.unique_pools = static_cast<int>(seen_pools.size());
    return out;
}

struct ExactPermutationMoments {
    bool ok = false;
    double mean = 0.0;
    double variance = 0.0;
    double skewness = 0.0;
    std::vector<double> centered_a;
};

ExactPermutationMoments exact_moments(const std::vector<double>& a,
                                      const std::vector<double>& g) {
    ExactPermutationMoments result;
    const int m = static_cast<int>(a.size());
    if (m < 6 || g.size() != a.size()) return result;

    const double mean_a = std::accumulate(a.begin(), a.end(), 0.0) / m;
    const double mean_g = std::accumulate(g.begin(), g.end(), 0.0) / m;
    result.centered_a.resize(a.size());
    double a2 = 0.0, a3 = 0.0, b2 = 0.0, b3 = 0.0;
    for (int i = 0; i < m; ++i) {
        const double ac = a[i] - mean_a;
        const double gc = g[i] - mean_g;
        result.centered_a[i] = ac;
        a2 += ac * ac;
        a3 += ac * ac * ac;
        b2 += gc * gc;
        b3 += gc * gc * gc;
    }
    const double variance = a2 * b2 / static_cast<double>(m - 1);
    if (!(variance > kTiny) || !std::isfinite(variance)) return result;

    const double third = a3 * b3 /
        (static_cast<double>(m - 1) * static_cast<double>(m - 2));
    result.mean = static_cast<double>(m) * mean_a * mean_g;
    result.variance = variance;
    result.skewness = third / std::pow(variance, 1.5);
    if (!std::isfinite(result.skewness)) result.skewness = 0.0;
    result.skewness = std::clamp(result.skewness, -2.5, 2.5);
    result.ok = true;
    return result;
}

double cornish_fisher_normal_score(double raw_z, double skewness) {
    // First-order Edgeworth probability match. For a variable with standardized
    // third cumulant gamma, F(z) ~= Phi(z) + gamma(1-z^2)phi(z)/6, hence the
    // Gaussian-equivalent quantile is z + gamma(1-z^2)/6.
    const double corrected = raw_z + skewness * (1.0 - raw_z * raw_z) / 6.0;
    return std::clamp(corrected, -12.0, 12.0);
}

std::vector<double> reference_for_band(
    const std::vector<std::vector<double>>& by_band,
    int target_band,
    int minimum_reference_size) {
    std::vector<double> reference;
    const int number_of_bands = static_cast<int>(by_band.size());
    for (int radius = 0; radius < number_of_bands; ++radius) {
        const int left = target_band - radius;
        const int right = target_band + radius;
        if (left >= 0) {
            reference.insert(reference.end(), by_band[left].begin(), by_band[left].end());
        }
        if (radius > 0 && right < number_of_bands) {
            reference.insert(reference.end(), by_band[right].begin(), by_band[right].end());
        }
        if (static_cast<int>(reference.size()) >= minimum_reference_size) break;
    }
    return reference;
}

}  // namespace

int permutation_size_band(int rows) {
    if (rows < 10) return 0;
    if (rows < 20) return 1;
    if (rows < 40) return 2;
    if (rows < 80) return 3;
    if (rows < 160) return 4;
    if (rows < 320) return 5;
    return 6;
}

PermutationStar build_permutation_star(Vertex& candidate,
                                       double null_height_units,
                                       int null_draws,
                                       std::uint64_t seed,
                                       bool allow_repeated_pool_rows) {
    PermutationStar result;
    const StarRows rows = collect_star_rows(candidate, null_height_units);
    result.rows = static_cast<int>(rows.a.size());
    result.unique_pools = rows.unique_pools;
    result.size_band = permutation_size_band(result.rows);
    if (result.rows < 6 || null_draws <= 0) return result;
    if (rows.repeated_pool_rows && !allow_repeated_pool_rows) return result;

    const ExactPermutationMoments moments = exact_moments(rows.a, rows.g);
    if (!moments.ok) return result;

    const double mean_g = std::accumulate(rows.g.begin(), rows.g.end(), 0.0) /
                          static_cast<double>(result.rows);
    double observed_centered = 0.0;
    for (int i = 0; i < result.rows; ++i) {
        observed_centered += moments.centered_a[i] * (rows.g[i] - mean_g);
    }
    result.raw_z = observed_centered / std::sqrt(moments.variance);
    result.skewness = moments.skewness;
    result.observed_z = cornish_fisher_normal_score(result.raw_z, result.skewness);

    std::vector<int> order(result.rows);
    std::iota(order.begin(), order.end(), 0);
    std::mt19937_64 rng(seed);
    result.null_z.reserve(null_draws);
    for (int draw = 0; draw < null_draws; ++draw) {
        std::shuffle(order.begin(), order.end(), rng);
        double numerator = 0.0;
        for (int i = 0; i < result.rows; ++i) {
            numerator += moments.centered_a[i] * (rows.g[order[i]] - mean_g);
        }
        const double raw = numerator / std::sqrt(moments.variance);
        result.null_z.push_back(cornish_fisher_normal_score(raw, result.skewness));
    }
    result.ok = std::isfinite(result.observed_z) && !result.null_z.empty();
    return result;
}

PooledPermutationResult pool_permutation_null(
    const std::vector<PermutationStar>& stars,
    int minimum_reference_size) {
    PooledPermutationResult result;
    result.p_values.assign(stars.size(), 1.0);
    result.gaussian_scores.assign(stars.size(), 0.0);
    result.reference_sizes.assign(stars.size(), 0);

    constexpr int kBands = 7;
    std::vector<std::vector<double>> by_band(kBands);
    for (const auto& star : stars) {
        if (!star.ok) continue;
        auto& band = by_band[std::clamp(star.size_band, 0, kBands - 1)];
        for (double value : star.null_z) {
            if (std::isfinite(value)) band.push_back(std::fabs(value));
        }
    }
    for (auto& band : by_band) std::sort(band.begin(), band.end());

    for (std::size_t i = 0; i < stars.size(); ++i) {
        const auto& star = stars[i];
        if (!star.ok) continue;
        std::vector<double> reference = reference_for_band(
            by_band, std::clamp(star.size_band, 0, kBands - 1),
            minimum_reference_size);
        if (reference.empty()) continue;
        std::sort(reference.begin(), reference.end());
        const double threshold = std::fabs(star.observed_z);
        const auto first = std::lower_bound(reference.begin(), reference.end(), threshold);
        const std::size_t exceedances = static_cast<std::size_t>(reference.end() - first);
        const double p = (static_cast<double>(exceedances) + 1.0) /
                         (static_cast<double>(reference.size()) + 1.0);
        const double clipped_p = std::clamp(p, 1e-15, 1.0);
        const double magnitude = inverse_standard_normal(1.0 - clipped_p / 2.0);
        result.p_values[i] = clipped_p;
        result.gaussian_scores[i] = star.observed_z < 0.0 ? -magnitude : magnitude;
        result.reference_sizes[i] = static_cast<int>(reference.size());
    }
    return result;
}

double inverse_standard_normal(double probability) {
    // Peter J. Acklam's rational approximation, with one Halley refinement.
    const double p = std::clamp(probability, 1e-15, 1.0 - 1e-15);
    constexpr double a1 = -3.969683028665376e+01;
    constexpr double a2 =  2.209460984245205e+02;
    constexpr double a3 = -2.759285104469687e+02;
    constexpr double a4 =  1.383577518672690e+02;
    constexpr double a5 = -3.066479806614716e+01;
    constexpr double a6 =  2.506628277459239e+00;
    constexpr double b1 = -5.447609879822406e+01;
    constexpr double b2 =  1.615858368580409e+02;
    constexpr double b3 = -1.556989798598866e+02;
    constexpr double b4 =  6.680131188771972e+01;
    constexpr double b5 = -1.328068155288572e+01;
    constexpr double c1 = -7.784894002430293e-03;
    constexpr double c2 = -3.223964580411365e-01;
    constexpr double c3 = -2.400758277161838e+00;
    constexpr double c4 = -2.549732539343734e+00;
    constexpr double c5 =  4.374664141464968e+00;
    constexpr double c6 =  2.938163982698783e+00;
    constexpr double d1 =  7.784695709041462e-03;
    constexpr double d2 =  3.224671290700398e-01;
    constexpr double d3 =  2.445134137142996e+00;
    constexpr double d4 =  3.754408661907416e+00;
    constexpr double plow = 0.02425;
    constexpr double phigh = 1.0 - plow;

    double x;
    if (p < plow) {
        const double q = std::sqrt(-2.0 * std::log(p));
        x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
            ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
    } else if (p > phigh) {
        const double q = std::sqrt(-2.0 * std::log(1.0-p));
        x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
    } else {
        const double q = p - 0.5;
        const double r = q*q;
        x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
            (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
    }
    const double cdf = 0.5 * std::erfc(-x / std::sqrt(2.0));
    const double pdf = std::exp(-0.5 * x * x) /
                       std::sqrt(2.0 * 3.14159265358979323846);
    x -= (cdf - p) / std::max(pdf, 1e-300);
    return x;
}
