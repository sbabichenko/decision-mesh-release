#pragma once

#include "dataset.h"
#include "run_config.h"
#include "surface.h"

int run_pipeline(const RunConfig& config, Dataset& dataset, const SurfaceSpec& surface);
