#include "dataset.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include "units.h"
#include <random>
#include <stdexcept>
#include <vector>

namespace {
constexpr double kBaseLogit = 0.0;
constexpr double kObservedCovariateCoefficient = 1.0;

Dataset load_real_data(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open data file: " + path);

    std::string line;
    if (!std::getline(file, line)) {
        throw std::runtime_error("data file is empty: " + path);
    }

    // DMESH_OFFSET5=1: csv carries a 5th column with a per-pool linear
    // predictor offset (known structure fitted outside the mesh). The mesh
    // then models only residual surface; gate and EB tax residual novelty.
    const bool use_off5 = std::getenv("DMESH_OFFSET5") != nullptr;
    const bool has_sigma2_pool = line.find("sigma2_pool") != std::string::npos;
    const bool has_n_eff = line.find("n_eff") != std::string::npos;
    std::vector<std::array<double, 7>> rows;
    double total_events = 0.0;
    double total_exposure = 0.0;
    double min_wala = 1e9, max_wala = -1e9;
    double min_wac = 1e9, max_wac = -1e9;
    int line_number = 1;
    while (std::getline(file, line)) {
        ++line_number;
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream parser(line);
        double wala = 0.0, wac = 0.0, exposure = 0.0, events = 0.0, off5 = 0.0;
        double sigma2_pool = -1.0, n_eff = -1.0;
        if (!(parser >> wala >> wac >> exposure >> events)) {
            throw std::runtime_error("invalid data row " + std::to_string(line_number) +
                                     " in " + path);
        }
        if (use_off5 && !(parser >> off5)) {
            throw std::runtime_error("DMESH_OFFSET5 set but row " +
                                     std::to_string(line_number) + " has no 5th column");
        } else if (!use_off5 && has_sigma2_pool && !(parser >> sigma2_pool)) {
            throw std::runtime_error("sigma2_pool header present but row " +
                                     std::to_string(line_number) + " has no value");
        }
        if (!use_off5 && has_n_eff && !(parser >> n_eff)) {
            throw std::runtime_error("n_eff header present but row " +
                                     std::to_string(line_number) + " has no value");
        }
        rows.push_back({wala, wac, exposure, events, off5, sigma2_pool, n_eff});
        total_events += events;
        total_exposure += exposure;
        min_wala = std::min(min_wala, wala);
        max_wala = std::max(max_wala, wala);
        min_wac = std::min(min_wac, wac);
        max_wac = std::max(max_wac, wac);
    }

    if (rows.empty() || total_exposure <= 0.0 || total_events <= 0.0 ||
        total_events >= total_exposure || max_wala <= min_wala || max_wac <= min_wac) {
        throw std::runtime_error("data file has invalid or degenerate contents: " + path);
    }

    const double pooled_probability = total_events / total_exposure;
    const double offset = std::log(pooled_probability / (1.0 - pooled_probability));

    Dataset dataset;
    dataset.observations.reserve(rows.size());
    for (const auto& row : rows) {
        const double x = (row[0] - min_wala) / (max_wala - min_wala);
        const double y = (row[1] - min_wac) / (max_wac - min_wac);
        const int n = static_cast<int>(row[2]);
        const int k = static_cast<int>(row[3]);
        double p0 = pooled_probability;
        double obs_offset = offset;
        if (use_off5) {
            obs_offset = offset + row[4];
            p0 = 1.0 / (1.0 + std::exp(-obs_offset));
        }
        const double information = n * p0 * (1.0 - p0);
        dataset.observations.push_back({x, y, information / kVarianceScale,
                                        from_log_odds((k - n * p0) / information),
                                        k, n, dataset.pool_count++, obs_offset, 0.0, 0.0,
                                        row[5], row[6]});
    }

    std::fprintf(stderr,
                 "[data] %zu pools, pooled rate %.4f, WALA [%g,%g], WAC [%g,%g]\n",
                 rows.size(), pooled_probability, min_wala, max_wala, min_wac, max_wac);
    dataset.truth_sigma_u = 0.0;
    return dataset;
}

Dataset simulate_data(const RunConfig& config, const SurfaceSpec& surface) {
    std::mt19937_64 location_rng(config.simulation.seed);
    std::mt19937_64 outcome_rng(config.simulation.noise_seed.value_or(config.simulation.seed));
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);

    Dataset dataset;
    dataset.truth_sigma_u = config.simulation.sigma_u;

    // Controlled benchmark: essentially uniform support over the full WALA-WAC
    // rectangle, with no attrition.  Every pool is observed at ages 1..60 and
    // WAC is drawn independently over its full range.  This is the default
    // synthetic experiment because it cleanly reveals whether the adaptive mesh
    // refines where the underlying function becomes more oscillatory.
    const int pools = static_cast<int>(std::round(936.0 * config.simulation.data_scale));
    for (int j = 0; j < pools; ++j, ++dataset.pool_count) {
        const double wac = 2.6 + 5.0 * uniform(location_rng);
        const int loans = static_cast<int>(std::exp(
            uniform(location_rng) * (std::log(3000.0) - std::log(50.0)) + std::log(50.0)));
        const double observed_covariate = normal(location_rng);
        const double persistent_effect = config.simulation.sigma_u
            * std::sqrt(config.simulation.permanent_variance_share) * normal(outcome_rng);

        for (int age = 1; age <= 60; ++age) {
            const double white_effect = config.simulation.sigma_u
                * std::sqrt(1.0 - config.simulation.permanent_variance_share) * normal(outcome_rng);
            const double pool_effect = persistent_effect + white_effect;

            // Mesh coordinates remain normalized on [0,1].  For interpretation,
            // the horizontal physical coordinate is x = WALA/12 in [1/12,5].
            const double x_raw = -4.0 + 8.0 * age / 60.0;
            const double y_raw = -4.0 + 8.0 * (wac - 2.6) / 5.0;
            const double systematic = surface.evaluate(x_raw, y_raw);
            const double true_logit = kBaseLogit + systematic
                + kObservedCovariateCoefficient * observed_covariate + pool_effect;
            const double true_probability = 1.0 / (1.0 + std::exp(-true_logit));
            std::binomial_distribution<int> binomial(loans, true_probability);
            const int events = binomial(outcome_rng);

            const double offset = kBaseLogit
                + kObservedCovariateCoefficient * observed_covariate;
            const double offset_probability = 1.0 / (1.0 + std::exp(-offset));
            const double information = loans * offset_probability * (1.0 - offset_probability);
            dataset.observations.push_back({(x_raw + 4.0) / 8.0,
                                            (y_raw + 4.0) / 8.0,
                                            information / kVarianceScale,
                                            (events - loans * offset_probability) / information * kHeightScale,
                                            events, loans, dataset.pool_count, offset,
                                            systematic * kHeightScale, pool_effect, -1.0, -1.0});
        }
    }
    return dataset;
}
}  // namespace

Dataset build_dataset(const RunConfig& config, const SurfaceSpec& surface) {
    if (config.simulation.data_file) return load_real_data(*config.simulation.data_file);
    return simulate_data(config, surface);
}
