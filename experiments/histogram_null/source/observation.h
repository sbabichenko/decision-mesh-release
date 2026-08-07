#pragma once

struct Observation {
    double x = 0.0;
    double y = 0.0;
    double working_weight = 0.0;
    double working_response = 0.0;
    int events = 0;
    int exposure = 0;
    int pool_id = 0;
    double offset = 0.0;
    double true_surface_scaled = 0.0;
    double latent_effect = 0.0;
    // Optional null-law metadata supplied by benchmark/production CSVs.
    // pool_variance is logit^2. effective_exposure is diagnostic only.
    double pool_variance = -1.0;
    double effective_exposure = -1.0;
};
