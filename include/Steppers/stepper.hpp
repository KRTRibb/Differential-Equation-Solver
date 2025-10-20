#pragma once
#include "types.hpp"

namespace diffeq {

class Stepper {
    public:
        virtual void step(const RHS& f, double& t, Vec& y, double h) = 0;

        virtual int GetOrder() const = 0;
        virtual const char* GetName() const = 0;

        virtual ~Stepper() = default;
};

}