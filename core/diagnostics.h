#pragma once
#include <vector>
struct DecisionMesh;
struct Observation;

[[nodiscard]] double pm_sigma_total(DecisionMesh* mesh,
    const std::vector<Observation>& observations,
    const std::vector<double>& scaled_x,
    const std::vector<double>& scaled_y);

[[nodiscard]] double leverage_correction(DecisionMesh* mesh,
    const std::vector<Observation>& observations,
    const std::vector<double>& scaled_x,
    const std::vector<double>& scaled_y,
    int pool_count,
    double naive_sigma,
    double& signal_denominator_ratio,
    double& noise_correction_ratio,
    int leverage_probes = 6);
