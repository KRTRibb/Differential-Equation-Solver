#include "Steppers/rk4_stepper.hpp"

using namespace diffeq;

void RK4Stepper::step(const RHS& f, double& t, Vec& y, double h) {
    Vec k1 = f(t, y);

    Vec y2(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y2[i] = y[i] + 0.5 * h * k1[i];
    }
    Vec k2 = f(t + 0.5 * h, y2);

    Vec y3(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y3[i] = y[i] + 0.5 * h * k2[i];
    }
    Vec k3 = f(t + 0.5 * h, y3);

    Vec y4(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y4[i] = y[i] + h * k3[i];
    }  
    Vec k4 = f(t + h, y4);

    for (size_t i = 0; i < y.size(); ++i) {
        y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0*k3[i] + k4[i]);
    }

    t += h;
}