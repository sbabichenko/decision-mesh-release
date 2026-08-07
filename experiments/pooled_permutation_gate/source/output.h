#pragma once

#include "observation.h"
#include "run_config.h"
#include "surface.h"

#include <vector>

class DecisionMesh;

void write_edf_dump(const OutputConfig& output,
                    const DecisionMesh& mesh,
                    int observation_count);

void write_figure_dumps(const OutputConfig& output,
                        const DecisionMesh& mesh,
                        const std::vector<Observation>& observations,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& original_x,
                        const std::vector<double>& original_y,
                        const std::vector<double>& x_map,
                        const std::vector<double>& y_map,
                        const SurfaceSpec& surface,
                        int pool_count,
                        double corrected_pool_variance);
