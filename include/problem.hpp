#pragma once
#include "types.hpp"
#include <utility>

namespace diffeq {

struct IVPProblem {
    RHS f;
    Vec y0;
    double t0;

    IVPProblem() = default;

    // Templated constructor to accept any callable
    template<typename F>
    IVPProblem(F&& func, const Vec& y_init, double t_start)
        : f(std::forward<F>(func)), y0(y_init), t0(t_start) {}
};

struct BVPProblem {
    RHS f;
    double t0;
    double t_end;
    Vec y0;
    Vec y_end;

    BVPProblem() = default;

    // Templated constructor for any callable
    template<typename F>
    BVPProblem(F&& func, double t_start, double t_finish,
               const Vec& y_init, const Vec& y_final)
        : f(std::forward<F>(func)), t0(t_start), t_end(t_finish),
          y0(y_init), y_end(y_final) {}
};

}
