#include "dataset.h"
#include "pipeline.h"
#include "run_config.h"
#include "surface.h"

#include <cstdio>
#include <exception>

int main(int argc, char** argv) {
    RunConfig config;
    try {
        config = parse_run_config(argc, argv);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "error: %s\n", error.what());
        print_usage(argv[0]);
        return 2;
    }

    const SurfaceSpec surface{config.simulation.chirp_surface, config.simulation.amplitude};
    try {
        Dataset dataset = build_dataset(config, surface);
        return run_pipeline(config, dataset, surface);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "error: %s\n", error.what());
        return 1;
    }
}
