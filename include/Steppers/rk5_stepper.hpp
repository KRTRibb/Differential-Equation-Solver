#pragma once
#include "stepper.hpp"

namespace diffeq {

class RK5Stepper : public Stepper {
    public:
        void step(const RHS& f, double& t, Vec& y, double h) override;

        int GetOrder() const override {return order;};
        const char* GetName() const override {return name;};

    private:
        const char* name = "Runge-Kutta 5";
        const int order = 5;
};

}