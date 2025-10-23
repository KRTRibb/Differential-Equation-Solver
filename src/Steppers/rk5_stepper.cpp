#include "Steppers/rk5_stepper.hpp"

using namespace diffeq;

void RK5Stepper::step(const RHS& f, double& t, Vec& y, double h) {
    Vec k1 = f(t, y);

    Vec y2(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        y2[i] = y[i] + 0.25 * h * k1[i];
    Vec k2 = f(t + 0.25 * h, y2);

    Vec y3(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        y3[i] = y[i] + (h / 8.0) * (k1[i] + k2[i]);
    Vec k3 = f(t + 0.25 * h, y3);

    Vec y4(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        y4[i] = y[i] + h * (-0.5 * k2[i] + k3[i]);
    Vec k4 = f(t + 0.5 * h, y4);

    Vec y5(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        y5[i] = y[i] + h * ((3.0/16.0) * k1[i] + (9.0/16.0) * k4[i]);
    Vec k5 = f(t + 0.75 * h, y5);

    Vec y6(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        y6[i] = y[i] + h * ((-3.0/7.0) * k1[i] + (2.0/7.0) * k2[i] +
                             (12.0/7.0) * k3[i] - (12.0/7.0) * k4[i] +
                             (8.0/7.0) * k5[i]);
    Vec k6 = f(t + h, y6);

    for (size_t i = 0; i < y.size(); ++i)
        y[i] += (h / 90.0) * (7.0*k1[i] + 32.0*k3[i] + 12.0*k4[i] + 32.0*k5[i] + 7.0*k6[i]);

    t += h;
}
