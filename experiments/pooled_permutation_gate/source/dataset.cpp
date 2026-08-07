#include "dataset.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>
#include <tuple>
#include <unordered_map>
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

struct SplitUnit {
    int pool_id = -1;
    double x = 0.0;
    double y = 0.0;
    double exposure = 0.0;
    std::uint64_t tie = 0;
};

void assign_split_mask(const RunConfig& config, Dataset& dataset) {
    const int n = static_cast<int>(dataset.observations.size());
    dataset.training_mask.assign(n, 1);
    if (!config.fit.heldout_split || n == 0) return;

    // External mask replay (DMESH_SPLIT_MASK): one 0/1 per data row, 1 = train.
    if (const char* mask_path = std::getenv("DMESH_SPLIT_MASK")) {
        std::ifstream mask_file(mask_path);
        if (!mask_file) throw std::runtime_error("DMESH_SPLIT_MASK: cannot open file");
        int value = 0, count = 0;
        while (mask_file >> value) {
            if (count >= n) throw std::runtime_error("DMESH_SPLIT_MASK: more entries than rows");
            dataset.training_mask[count++] = (value != 0);
        }
        if (count != n) throw std::runtime_error("DMESH_SPLIT_MASK: entry count != row count");
        return;
    }

    if (config.fit.split_mode == "parity") {
        for (int i = 0; i < n; ++i) dataset.training_mask[i] = (i % 2 == 0);
    } else {
        // Split whole pools, not rows. This avoids leakage on longitudinal
        // synthetic data while reducing to one unit per row for the Ginnie
        // snapshot. A median k-d partition forms small adaptive spatial blocks;
        // exposure-ranked pairs inside each block are sent to opposite sides.
        const int pool_count = std::max(dataset.pool_count, 0);
        std::vector<double> sx(pool_count, 0.0), sy(pool_count, 0.0), sw(pool_count, 0.0);
        std::vector<std::vector<int>> rows(pool_count);
        for (int i = 0; i < n; ++i) {
            const Observation& o = dataset.observations[i];
            if (o.pool_id < 0 || o.pool_id >= pool_count) continue;
            const double w = std::max(o.exposure, 1);
            sx[o.pool_id] += w * o.x;
            sy[o.pool_id] += w * o.y;
            sw[o.pool_id] += w;
            rows[o.pool_id].push_back(i);
        }

        std::mt19937_64 rng(config.simulation.seed);
        std::vector<SplitUnit> units;
        units.reserve(pool_count);
        for (int pool = 0; pool < pool_count; ++pool) {
            if (rows[pool].empty()) continue;
            units.push_back({pool, sx[pool] / sw[pool], sy[pool] / sw[pool], sw[pool], rng()});
        }
        std::vector<int> order(units.size());
        std::iota(order.begin(), order.end(), 0);
        std::vector<unsigned char> pool_train(pool_count, 0);
        std::size_t block_count = 0;
        double max_block_exposure_imbalance = 0.0;

        std::function<void(int,int)> split_block = [&](int begin, int end) {
            const int count = end - begin;
            if (count <= config.fit.split_block_size) {
                ++block_count;
                std::sort(order.begin() + begin, order.begin() + end,
                          [&](int a, int b) {
                              if (units[a].exposure != units[b].exposure)
                                  return units[a].exposure > units[b].exposure;
                              return units[a].tie < units[b].tie;
                          });
                // Greedy equal-cardinality partition minimizes exposure
                // imbalance inside the local block. A random block-wide label
                // flip supplies seeded variation without sacrificing balance.
                const int target_a = (count + 1) / 2;
                const int target_b = count / 2;
                std::vector<unsigned char> side_a(count, 0);
                double exp_a = 0.0, exp_b = 0.0;
                int na = 0, nb = 0;
                for (int j = 0; j < count; ++j) {
                    const int u = order[begin + j];
                    bool to_a;
                    if (na >= target_a) to_a = false;
                    else if (nb >= target_b) to_a = true;
                    else if (std::fabs(exp_a - exp_b) <= 1e-12) to_a = (rng() & 1ULL) != 0;
                    else to_a = exp_a < exp_b;
                    side_a[j] = to_a ? 1 : 0;
                    if (to_a) { exp_a += units[u].exposure; ++na; }
                    else { exp_b += units[u].exposure; ++nb; }
                }
                const bool flip_labels = (rng() & 1ULL) != 0;
                for (int j = 0; j < count; ++j) {
                    const int u = order[begin + j];
                    const bool train = (side_a[j] != 0) ^ flip_labels;
                    pool_train[units[u].pool_id] = train ? 1 : 0;
                }
                const double train_exp = flip_labels ? exp_b : exp_a;
                const double held_exp = flip_labels ? exp_a : exp_b;
                const double total = train_exp + held_exp;
                if (total > 0.0)
                    max_block_exposure_imbalance = std::max(
                        max_block_exposure_imbalance,
                        std::fabs(train_exp - held_exp) / total);
                return;
            }

            double xmin = 1.0, xmax = 0.0, ymin = 1.0, ymax = 0.0;
            for (int j = begin; j < end; ++j) {
                const SplitUnit& u = units[order[j]];
                xmin = std::min(xmin, u.x); xmax = std::max(xmax, u.x);
                ymin = std::min(ymin, u.y); ymax = std::max(ymax, u.y);
            }
            const bool x_axis = (xmax - xmin) >= (ymax - ymin);
            const int mid = begin + count / 2;
            std::nth_element(order.begin() + begin, order.begin() + mid, order.begin() + end,
                             [&](int a, int b) {
                                 const double ca = x_axis ? units[a].x : units[a].y;
                                 const double cb = x_axis ? units[b].x : units[b].y;
                                 if (ca != cb) return ca < cb;
                                 const double oa = x_axis ? units[a].y : units[a].x;
                                 const double ob = x_axis ? units[b].y : units[b].x;
                                 if (oa != ob) return oa < ob;
                                 return units[a].tie < units[b].tie;
                             });
            split_block(begin, mid);
            split_block(mid, end);
        };
        split_block(0, static_cast<int>(order.size()));

        for (int pool = 0; pool < pool_count; ++pool)
            for (int i : rows[pool]) dataset.training_mask[i] = pool_train[pool];

        std::fprintf(stderr,
                     "[split] local seed %u, block<=%d, %zu blocks, max block exposure imbalance %.3f\n",
                     config.simulation.seed, config.fit.split_block_size, block_count,
                     max_block_exposure_imbalance);
    }

    long train_rows = 0, held_rows = 0;
    double train_n = 0.0, held_n = 0.0, train_k = 0.0, held_k = 0.0;
    for (int i = 0; i < n; ++i) {
        const Observation& o = dataset.observations[i];
        if (dataset.training_mask[i]) {
            ++train_rows; train_n += o.exposure; train_k += o.events;
        } else {
            ++held_rows; held_n += o.exposure; held_k += o.events;
        }
    }
    std::fprintf(stderr,
                 "[split] mode %s: train %ld rows / %.0f exposure / rate %.6f; "
                 "held %ld rows / %.0f exposure / rate %.6f\n",
                 config.fit.split_mode.c_str(), train_rows, train_n, train_k / std::max(train_n, 1.0),
                 held_rows, held_n, held_k / std::max(held_n, 1.0));
}
}  // namespace

Dataset build_dataset(const RunConfig& config, const SurfaceSpec& surface) {
    Dataset dataset = config.simulation.data_file
        ? load_real_data(*config.simulation.data_file)
        : simulate_data(config, surface);
    assign_split_mask(config, dataset);
    return dataset;
}
