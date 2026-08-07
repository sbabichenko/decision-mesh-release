#pragma once
#include <vector>

struct LfdrResult {
    std::vector<double> lfdr;
    double empirical_mean = 0.0;
    double empirical_sd = 1.0;
    double pi0 = 1.0;
    bool used_lindsey = false;
};

[[nodiscard]] LfdrResult lindsey_lfdr(std::vector<double> z,
                                        bool free_empirical_null = false);
[[nodiscard]] double normal_survival(double z);
