#pragma once

#include <cmath>

struct SurfaceSpec {
    bool chirp = false;
    double amplitude = 2.0;

    [[nodiscard]] double evaluate(double x, double y) const {
        // x,y are the simulator's internal coordinates on [-4,4].
        // Convert x to z in [0,1], then accelerate its phase smoothly so
        // horizontal frequency increases toward the right side of the domain.
        const double z = (x + 4.0) / 8.0;
        const double z_warp = (0.2 + 0.8 * z) * z;
        const double x_eff = 5.0 * z_warp;
        const double xr = -4.0 + 8.0 * x_eff;

        const double coarse =
            std::exp(0.80 * std::cos(0.72 * xr + 0.44 * y + 0.25))
            - 1.1665149229;
        const double fine =
            std::exp(0.58 * std::cos(2.65 * (-0.36 * xr + 0.93 * y) - 0.40))
            - 1.0859474;
        return amplitude * (0.62 * coarse + 0.38 * fine);
    }
};
