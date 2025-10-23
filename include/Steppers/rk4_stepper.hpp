#pragma once
#include "stepper.hpp"

namespace diffeq {


class RK4Stepper : public Stepper {
    public:
        void step(const RHS& f, double& t, Vec& y, double h) override;

        int GetOrder() const override {return order;};
        const char* GetName() const override {return name;};

    private:
        const char* name = "Runge-Kutta 4";
        const int order = 3;
};

}