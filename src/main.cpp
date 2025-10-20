#include <iostream>
#include <cmath>
#include "solver.hpp"
#include "Steppers/euler_stepper.hpp"
#include "errorfunc/funcs.hpp"
#include "errorfunc/metric.hpp"

using namespace diffeq;
using namespace diffeq::errorfunc;

int main() {
    // Choose a stepper (can swap easily)
    EulerStepper euler;
    Solver solver(euler);

    // ---- TEST 1 ----
    // dy/dt = -y, y(0) = 1 -> exact: y(t) = exp(-t)
    IVPProblem prob1;
    prob1.f = [](double t, const Vec& y) {
        Vec dydt;
        dydt.resize(1);
        dydt[0] = -y[0];
        return dydt;
    };
    prob1.y0 = {1.0};
    prob1.t0 = 0.0;

    auto exact1 = [](double t) { return Vec{std::exp(-t)}; };

    std::cout << "\n--- TEST 1: dy/dt = -y ---\n";
    solver.demoRun(prob1, 1.0, 0.1, exact1, RMS);
    solver.convergenceTestDemoRun(prob1, 1.0, exact1, 1, 10, RMS);

    // ---- TEST 2 ----
    // dy/dt = cos(t), y(0) = 0 -> exact: y(t) = sin(t)
    IVPProblem prob2;
    prob2.f = [](double t, const Vec& y) {
        Vec dydt;
        dydt.resize(1);
        dydt[0] = std::cos(t);
        return dydt;
    };
    prob2.y0 = {0.0};
    prob2.t0 = 0.0;

    auto exact2 = [](double t) { return Vec{std::sin(t)}; };

    std::cout << "\n--- TEST 2: dy/dt = cos(t) ---\n";
    solver.demoRun(prob2, M_PI, 0.1, exact2, RMS);
    solver.convergenceTestDemoRun(prob2, M_PI, exact2, 1, 10, RMS);

    // ---- TEST 3 ----
    // System dy/dt = A*y where A = [[0,1],[-1,0]] (circular motion)
    // Exact solution: y1 = cos(t), y2 = -sin(t)
    IVPProblem prob3;
    prob3.f = [](double t, const Vec& y) {
        Vec dydt;
        dydt.resize(2);
        dydt[0] = y[1];
        dydt[1] = -y[0];

        return dydt;
    };
    prob3.y0 = {1.0, 0.0};
    prob3.t0 = 0.0;

    double t_end = 2 * M_PI;
    double h = 0.1;

    auto exact3 = [](double t) { return Vec{std::cos(t), -std::sin(t)}; };

    std::cout << "\n--- TEST 3: Circular Motion System ---\n";
    solver.demoRun(prob3, t_end, h, exact3, REL_L2);
    solver.convergenceTestDemoRun(prob3, t_end, exact3, 1, 10, REL_L2);

    std::cout << "\nAll tests complete.\n";
    return 0;
}
