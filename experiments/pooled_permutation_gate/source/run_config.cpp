#include "run_config.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

bool env_flag(const char* name) {
    const char* value = std::getenv(name);
    if (value == nullptr) return false;
    const std::string text(value);
    return text.empty() || (text != "0" && text != "false" && text != "FALSE");
}

std::optional<std::string> env_string(const char* name) {
    const char* value = std::getenv(name);
    if (value == nullptr || *value == '\0') return std::nullopt;
    return std::string(value);
}


void reject_retired_option(const char* name, const char* replacement) {
    if (std::getenv(name) == nullptr) return;
    std::string message = std::string(name) + " was removed";
    if (replacement && *replacement) message += std::string("; ") + replacement;
    throw std::runtime_error(message);
}

template <typename T>
T parse_number(const char* text, const char* label);

template <>
double parse_number<double>(const char* text, const char* label) {
    char* end = nullptr;
    errno = 0;
    const double value = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !std::isfinite(value))
        throw std::runtime_error(std::string("invalid ") + label + ": " + text);
    return value;
}

template <>
int parse_number<int>(const char* text, const char* label) {
    char* end = nullptr;
    errno = 0;
    const long value = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
        throw std::runtime_error(std::string("invalid ") + label + ": " + text);
    return static_cast<int>(value);
}

template <typename T>
T env_number(const char* name, T fallback) {
    const char* value = std::getenv(name);
    return value == nullptr ? fallback : parse_number<T>(value, name);
}

}  // namespace

RunConfig parse_run_config(int argc, char** argv) {
    reject_retired_option("DMESH_LEGACY_EB", "the zero-mean surplus prior is now the only implementation");
    reject_retired_option("DMESH_TAU_ALL_ACTIVE", "tau is now always estimated from gate admissions");
    reject_retired_option("DMESH_PM_LEGACY", "the deterministic fixed-point PM solver is now the only implementation");

    RunConfig config;
    if (argc > 1) config.simulation.sigma_u = parse_number<double>(argv[1], "sigma_u");
    if (argc > 2) config.fit.fdr_level = parse_number<double>(argv[2], "q");
    if (argc > 3) config.simulation.data_scale = parse_number<int>(argv[3], "scale");
    if (argc > 4) config.fit.points_per_face = parse_number<int>(argv[4], "points_per_face");
    if (argc > 5) {
        const int parsed_seed = parse_number<int>(argv[5], "seed");
        if (parsed_seed < 0) throw std::runtime_error("seed must be nonnegative");
        config.simulation.seed = static_cast<unsigned>(parsed_seed);
    }
    if (argc > 6) throw std::runtime_error("too many command-line arguments");

    config.simulation.chirp_surface = env_flag("DMESH_CHIRP");
    config.simulation.amplitude = env_number<double>("DMESH_AMP", config.simulation.amplitude);
    config.simulation.permanent_variance_share = env_number<double>("DMESH_PERM", config.simulation.permanent_variance_share);
    config.simulation.attrition_rate = env_number<double>("DMESH_ATTR", config.simulation.attrition_rate);
    if (const auto seed = env_string("DMESH_NOISESEED")) {
        const int parsed = parse_number<int>(seed->c_str(), "DMESH_NOISESEED");
        if (parsed < 0) throw std::runtime_error("DMESH_NOISESEED must be nonnegative");
        config.simulation.noise_seed = static_cast<unsigned>(parsed);
    }
    config.simulation.data_file = env_string("DMESH_DATA");

    config.fit.learned_metric = env_flag("DMESH_METRIC");
    config.fit.oracle_vertex_budget = env_number<int>("DMESH_ORACLE", 0);
    config.fit.disable_eb = env_flag("DMESH_NOEB");
    config.fit.heldout_split = env_flag("DMESH_SPLIT");
    config.fit.split_mode = env_string("DMESH_SPLIT_MODE").value_or(config.fit.split_mode);
    config.fit.split_block_size = env_number<int>("DMESH_SPLIT_BLOCK_SIZE", config.fit.split_block_size);
    config.fit.parent_variance_scale = env_number<double>(
        "DMESH_PARENT_VAR_SCALE", config.fit.parent_variance_scale);
    config.diagnostics.dump_score_diagnostics = env_flag("DMESH_ZDIAG");
    config.diagnostics.lag_match_bands = env_number<int>("DMESH_BANDS", config.diagnostics.lag_match_bands);
    config.diagnostics.max_lag = env_number<int>("DMESH_MAXLAG", config.diagnostics.max_lag);
    config.diagnostics.oracle_shared = env_flag("DMESH_ORACLE_SHARED");
    config.diagnostics.run_height_test = env_flag("DMESH_HTEST");
    config.output.edf_dump_file = env_string("DMESH_EDFDUMP");
    config.output.prefix = env_string("DMESH_DUMP");
    config.diagnostics.leverage_probes = env_number<int>("DMESH_LEVERAGE_PROBES", config.diagnostics.leverage_probes);
    config.diagnostics.max_gate_rounds = env_number<int>("DMESH_MAX_GATE_ROUNDS", config.diagnostics.max_gate_rounds);
    config.diagnostics.stage_count = env_number<int>("DMESH_STAGE_COUNT", config.diagnostics.stage_count);
    config.diagnostics.stop_after_adaptation = env_flag("DMESH_STOP_AFTER_ADAPTATION");

    if (!(config.simulation.sigma_u >= 0.0)) throw std::runtime_error("sigma_u must be nonnegative");
    if (!(config.fit.fdr_level > 0.0 && config.fit.fdr_level < 1.0))
        throw std::runtime_error("q must lie strictly between 0 and 1");
    if (config.simulation.data_scale <= 0) throw std::runtime_error("scale must be positive");
    if (config.fit.points_per_face <= 0) throw std::runtime_error("points_per_face must be positive");
    if (config.fit.split_mode != "local" && config.fit.split_mode != "parity")
        throw std::runtime_error("DMESH_SPLIT_MODE must be local or parity");
    if (config.fit.split_block_size < 4 || config.fit.split_block_size > 256)
        throw std::runtime_error("DMESH_SPLIT_BLOCK_SIZE must lie in [4,256]");
    if (!(config.fit.parent_variance_scale >= 0.0 &&
          std::isfinite(config.fit.parent_variance_scale)))
        throw std::runtime_error("DMESH_PARENT_VAR_SCALE must be finite and nonnegative");
    if (!(config.simulation.permanent_variance_share >= 0.0 && config.simulation.permanent_variance_share <= 1.0))
        throw std::runtime_error("DMESH_PERM must lie in [0,1]");
    if (!(config.simulation.attrition_rate >= 0.0)) throw std::runtime_error("DMESH_ATTR must be nonnegative");
    if (config.fit.oracle_vertex_budget < 0) throw std::runtime_error("DMESH_ORACLE must be nonnegative");
    if (config.diagnostics.lag_match_bands <= 0) throw std::runtime_error("DMESH_BANDS must be positive");
    if (config.diagnostics.max_lag <= 0) throw std::runtime_error("DMESH_MAXLAG must be positive");
    config.diagnostics.leverage_probes = std::clamp(config.diagnostics.leverage_probes, 1, 64);
    if (config.diagnostics.max_gate_rounds <= 0 || config.diagnostics.max_gate_rounds > 80)
        throw std::runtime_error("DMESH_MAX_GATE_ROUNDS must lie in [1,80]");
    if (config.diagnostics.stage_count <= 0 || config.diagnostics.stage_count > 3)
        throw std::runtime_error("DMESH_STAGE_COUNT must lie in [1,3]");

    return config;
}

void print_usage(const char* program_name) {
    std::fprintf(stderr,
                 "usage: %s [sigma_u=0.30] [q=0.10] [scale=1] "
                 "[points_per_face=12] [seed=7]\n",
                 program_name);
}
