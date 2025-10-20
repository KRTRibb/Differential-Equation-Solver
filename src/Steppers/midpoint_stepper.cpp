#include "Steppers/midpoint_stepper.hpp"

using namespace diffeq;

void MidpointStepper::step(const RHS& f, double& t, Vec& y, double h) {
    Vec k1 = f(t, y);

    // compute y at midpoint
    Vec y_mid(y.size());
    for (size_t i = 0; i < y.size(); ++i){
        y_mid[i] = y[i] + 0.5 * h * k1[i];
    }

    Vec k2 = f(t + 0.5 * h, y_mid);

    for (size_t i = 0; i < y.size(); ++i) {
        y[i] += h * k2[i];
    }

    t += h;
}