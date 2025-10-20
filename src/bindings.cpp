#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include "solver.hpp"
#include "Steppers/stepper.hpp"
#include "problem.hpp"
#include "Steppers/euler_stepper.hpp"
#include "Steppers/midpoint_stepper.hpp"
#include "Steppers/rk3_stepper.hpp"
#include "errorfunc/metric.hpp"
#include "errorfunc/funcs.hpp"

namespace py = pybind11;
using namespace diffeq;

// Make std::vector<double> opaque for custom bindings
PYBIND11_MAKE_OPAQUE(std::vector<double>);

class PyStepper : public Stepper {
    using Stepper::Stepper;
    void step(const RHS& f, double& t, Vec& y, double h) override {
        PYBIND11_OVERRIDE_PURE(void, Stepper, step, f, t, y, h);
    }
    int GetOrder() const override {
        PYBIND11_OVERRIDE_PURE(int, Stepper, GetOrder, );
    }
    const char* GetName() const override {
        PYBIND11_OVERRIDE_PURE(const char*, Stepper, GetName, );
    }
};

PYBIND11_MODULE(diffeqpy, m) {
    m.doc() = "Python bindings for diffeq solver library";

    // ---- Vec bindings ----
    py::bind_vector<std::vector<double>>(m, "Vec")
        .def("__getitem__", [](const Vec &v, size_t i) { return v[i]; })
        .def("__setitem__", [](Vec &v, size_t i, double val) { v[i] = val; })
        .def("__len__", [](const Vec &v) { return v.size(); })
        .def(py::init([](const py::list &lst) {
            Vec v;
            v.reserve(py::len(lst));
            for (auto item : lst)
                v.push_back(py::cast<double>(item));
            return v;
        }));

    py::implicitly_convertible<py::list, std::vector<double>>();

    // ---- IVPProblem ----
    // Now f is expected to return a Vec
    py::class_<IVPProblem>(m, "IVPProblem")
        .def(py::init<const std::function<Vec(double, const Vec&)>&, const Vec&, double>(),
             py::arg("f"), py::arg("y0"), py::arg("t0"))
        .def(py::init([](const std::function<Vec(double, const Vec&)>& f,
                         const py::list& y0_list, double t0) {
            Vec y0;
            y0.reserve(py::len(y0_list));
            for (auto item : y0_list)
                y0.push_back(py::cast<double>(item));
            return IVPProblem{f, y0, t0};
        }), py::arg("f"), py::arg("y0"), py::arg("t0"))
        .def(py::init([](const std::function<Vec(double, const Vec&)>& f,
                         const py::array_t<double>& y0_array, double t0) {
            auto buf = y0_array.request();
            Vec y0(buf.size);
            std::memcpy(y0.data(), buf.ptr, buf.size * sizeof(double));
            return IVPProblem{f, y0, t0};
        }), py::arg("f"), py::arg("y0"), py::arg("t0"))
        .def_readwrite("f", &IVPProblem::f)
        .def_readwrite("y0", &IVPProblem::y0)
        .def_readwrite("t0", &IVPProblem::t0);

    // ---- ErrorMetric ----
    py::class_<errorfunc::ErrorMetric>(m, "ErrorMetric")
        .def(py::init<std::string, std::function<double(const Vec&, const Vec&)>>())
        .def_readwrite("name", &errorfunc::ErrorMetric::name)
        .def_readwrite("func", &errorfunc::ErrorMetric::func);
    m.attr("L2") = errorfunc::L2;
    m.attr("rms") = errorfunc::RMS;
    m.attr("absolute") = errorfunc::ABS;
    m.attr("relative") = errorfunc::REL;
    m.attr("max_norm") = errorfunc::MAX;
    m.attr("rel_L2") = errorfunc::REL_L2;


    // ---- IntegrationResult ----
    py::class_<IntegrationResult>(m, "IntegrationResult")
        .def_readonly("T", &IntegrationResult::T)
        .def_readonly("Y", &IntegrationResult::Y)
        .def_property_readonly("T_np", [](const IntegrationResult& self) {
            return py::array_t<double>(self.T.size(), self.T.data());
        })
        .def_property_readonly("Y_np", [](const IntegrationResult& self) {
            if (self.Y.empty()) return py::array_t<double>();
            size_t n_steps = self.Y.size();
            size_t n_dim = self.Y[0].size();
            auto arr = py::array_t<double>({n_steps, n_dim});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < n_steps; ++i)
                for (size_t j = 0; j < n_dim; ++j)
                    buf(i, j) = self.Y[i][j];
            return arr;
        })
        .def_readonly("exactY", &IntegrationResult::exactY)
        .def_readonly("errors", &IntegrationResult::errors)
        .def_readonly("error_func_name", &IntegrationResult::error_func_name)
        .def_readonly("n_steps", &IntegrationResult::n_steps)
        .def_readonly("h_used", &IntegrationResult::h_used)
        .def_readonly("total_time", &IntegrationResult::total_time)
        .def_readonly("final_error", &IntegrationResult::final_error);


    py::class_<ConvergenceTestResult>(m, "ConvergeceTestResult")
        .def_readonly("t_end", &ConvergenceTestResult::t_end)
        .def_readonly("h_vals", &ConvergenceTestResult::h_vals)
        .def_readonly("final_error", &ConvergenceTestResult::final_error)
        .def_readonly("final_error_ratios", &ConvergenceTestResult::final_error_ratios)
        .def_readonly("p_estimations", &ConvergenceTestResult::p_estimations)
        .def_readonly("error_name", &ConvergenceTestResult::error_name)
        .def_readonly("min_pow", &ConvergenceTestResult::min_pow)
        .def_readonly("max_pow", &ConvergenceTestResult::max_pow);

    // ---- Stepper ----
    py::class_<Stepper, PyStepper>(m, "Stepper")
        .def(py::init<>())
        .def("step", &Stepper::step)
        .def("GetOrder", &Stepper::GetOrder)
        .def("GetName", &Stepper::GetName);

    py::class_<EulerStepper, Stepper>(m, "EulerStepper")
        .def(py::init<>())
        .def("GetOrder", &EulerStepper::GetOrder)
        .def("GetName", &EulerStepper::GetName);

    py::class_<MidpointStepper, Stepper>(m, "MidpointStepper")
        .def(py::init<>())
        .def("GetOrder", &MidpointStepper::GetOrder)
        .def("GetName", &MidpointStepper::GetName);

    py::class_<RK3Stepper, Stepper>(m, "RK3Stepper")
        .def(py::init<>())
        .def("GetOrder", & RK3Stepper::GetOrder)
        .def("GetName", &RK3Stepper::GetName);

    // ---- Solver ----
    py::class_<Solver>(m, "Solver")
        .def(py::init<Stepper&>(), py::keep_alive<1, 2>())
        .def("integrateFixedSteps",
             py::overload_cast<const IVPProblem&, double, double>(&Solver::integrateFixedSteps, py::const_),
             py::arg("prob"), py::arg("t_end"), py::arg("h"))
        .def("integrateFixedSteps",
             py::overload_cast<const IVPProblem&, double, double, const std::function<Vec(double)>&, const errorfunc::ErrorMetric&>(&Solver::integrateFixedSteps, py::const_),
             py::arg("prob"), py::arg("t_end"), py::arg("h"), py::arg("exactSolution"), py::arg("metric") = errorfunc::L2)
        .def("printResults", &Solver::printResults)
        .def("GetApproximations", &Solver::GetApproximations)
        .def("GetExactVals", &Solver::GetExactVals)
        .def("GetErrors", &Solver::GetErrors)
        .def("convergenceTest", &Solver::convergenceTest)
        .def("printConvergenceTest", &Solver::printConvergenceTest);
}
