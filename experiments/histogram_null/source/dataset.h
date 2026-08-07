#pragma once

#include "observation.h"
#include "run_config.h"
#include "surface.h"

#include <vector>

struct Dataset {
    std::vector<Observation> observations;
    int pool_count = 0;
    double truth_sigma_u = 0.0;
};

Dataset build_dataset(const RunConfig& config, const SurfaceSpec& surface);
