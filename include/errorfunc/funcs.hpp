// Header file defining various error metrics for evaluating numerical methods
#pragma once

#include <vector>
#include <cmath>
#include <functional>
#include "types.hpp"

namespace diffeq::errorfunc {

// Computes the mean absolute error between two vectors.
// Formula: (1/n) * Σ |y_true[i] - y_approx[i]|
inline double absolute(const Vec& y_true, const Vec& y_approx) {
    double err = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i)
        err += std::abs(y_true[i] - y_approx[i]);
    return err / y_true.size();
}

// Computes the root mean square (RMS) error between two vectors.
// Formula: sqrt( (1/n) * Σ (y_true[i] - y_approx[i])² )
inline double rms(const Vec& y_true, const Vec& y_approx) {
    double err = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i)
        err += std::pow(y_true[i] - y_approx[i], 2);
    return std::sqrt(err / y_true.size());
}

// Computes the mean relative (fractional) error between two vectors.
// Formula: (1/n) * Σ |(y_true[i] - y_approx[i]) / (y_true[i] + ε)|
inline double relative(const Vec& y_true, const Vec& y_approx) {
    double err = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i)
        err += std::abs((y_true[i] - y_approx[i]) / (y_true[i] + 1e-12));
    return err / y_true.size();
}

// Computes the L2 norm (Euclidean distance) between two vectors.
// Formula: sqrt( Σ (y_true[i] - y_approx[i])² )
inline double l2norm(const Vec& y_true, const Vec& y_approx) {
    double sum = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i)
        sum += (y_true[i] - y_approx[i]) * (y_true[i] - y_approx[i]);
    return std::sqrt(sum);
}

// Computes the maximum absolute error across all components (L∞ norm).
// Formula: max_i |y_true[i] - y_approx[i]|
inline double maxnorm(const Vec& y_true, const Vec& y_approx) {
    double maxerr = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i)
        maxerr = std::max(maxerr, std::abs(y_true[i] - y_approx[i]));
    return maxerr;
}

// Computes the relative L2 norm between two vectors.
// Formula: sqrt( Σ (Δy[i])² / Σ (y_approx[i])² )
// Measures Euclidean distance normalized by the magnitude of y_approx.
inline double rel_l2norm(const Vec& y_true, const Vec& y_approx) {
    double num = 0.0, denom = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i) {
        num += (y_true[i] - y_approx[i]) * (y_true[i] - y_approx[i]);
        denom += y_approx[i] * y_approx[i];
    }
    return std::sqrt(num / denom);
}

} // namespace diffeq::errorfunc
