#pragma once
#include <vector>
#include <functional>

namespace diffeq {

    using Vec = std::vector<double>;
    using RHS = std::function<Vec(double t, const Vec& y)>;

    /*
    Later can add:

    Matrix alias (for Jacobians)

    Real type template (float/double/long double)

    Parallel vector implementations (Eigen, custom allocators, etc.)
    */
}