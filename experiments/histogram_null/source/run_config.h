#pragma once

#include <optional>
#include <string>

struct SimulationConfig {
    double sigma_u = 0.30;
    int data_scale = 1;
    unsigned seed = 7;
    bool chirp_surface = false;
    double amplitude = 2.0;
    double permanent_variance_share = 1.0;
    double attrition_rate = 0.05;
    std::optional<unsigned> noise_seed;
    std::optional<std::string> data_file;
};

struct FitConfig {
    double fdr_level = 0.10;
    int points_per_face = 12;
    bool learned_metric = false;
    int oracle_vertex_budget = 0;
    bool disable_eb = false;
    bool heldout_split = false;

    [[nodiscard]] bool oracle_mode() const { return oracle_vertex_budget > 0; }
    [[nodiscard]] bool use_eb() const { return !disable_eb && !oracle_mode(); }
};

struct DiagnosticConfig {
    bool dump_score_diagnostics = false;
    int lag_match_bands = 24;
    int max_lag = 24;
    bool oracle_shared = false;
    bool run_height_test = false;
    int leverage_probes = 6;
};

struct OutputConfig {
    std::optional<std::string> edf_dump_file;
    std::optional<std::string> prefix;
};

struct RunConfig {
    SimulationConfig simulation;
    FitConfig fit;
    DiagnosticConfig diagnostics;
    OutputConfig output;
};

RunConfig parse_run_config(int argc, char** argv);
void print_usage(const char* program_name);
