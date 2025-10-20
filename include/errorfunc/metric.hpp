#pragma once

#include "funcs.hpp"

namespace diffeq::errorfunc {

// Defining an errormetric struct that contains an error func and its associated name.
struct ErrorMetric {
    std::string name;
    std::function<double(const Vec&, const Vec&)> func;

    ErrorMetric(std::string n, std::function<double(const Vec&, const Vec&)> f)
        : name(std::move(n)), func(std::move(f)) {}
};

// predefined errors
inline const ErrorMetric L2{"L2 Norm error", l2norm};
inline const ErrorMetric RMS{"RMS Error", rms};
inline const ErrorMetric REL{"Mean Relative Error", relative};
inline const ErrorMetric ABS{"Mean Absolute Error", absolute};
inline const ErrorMetric MAX{"Max Norm Error", maxnorm};
inline const ErrorMetric REL_L2{"Relative L2 Norm Error", rel_l2norm};
}