#pragma once

#include "observation.h"
#include "run_config.h"
#include "surface.h"

#include <vector>

struct Dataset {
    std::vector<Observation> observations;
    int pool_count = 0;
    double truth_sigma_u = 0.0;
    // 1 for fitting rows, 0 for held-out rows. Empty only before split
    // assignment; run_pipeline always receives a mask of observations.size().
    std::vector<unsigned char> training_mask;
};

Dataset build_dataset(const RunConfig& config, const SurfaceSpec& surface);
