#pragma once
#include "types.hpp"

namespace diffeq {

// Inital value problem
struct IVPProblem {
RHS f;
Vec y0;
double t0;

IVPProblem() = default;

IVPProblem(const RHS& func, const Vec& y_init, double t_start)
        : f(func), y0(y_init), t0(t_start) {}
};

// Boundary value problem
struct BVPProblem {
RHS f;
double t0, t_end;
Vec y0;
Vec y_end;
};

}

