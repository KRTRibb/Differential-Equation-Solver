#include "pde.hpp"
#include <cassert>

using namespace diffeq;

RHS pde::make_heat_rhs(double alpha,
                        std::size_t n_interior,
                        double dx,
                        pde::BCType left_bc,
                        pde::BCType right_bc,
                        pde::BCFunc left_val,
                        pde::BCFunc right_val) {
    return [alpha, n_interior, dx, left_bc, right_bc, left_val = std::move(left_val), right_val = std::move(right_val)](double t, const Vec& y) -> Vec {
        assert(y.size() == n_interior);
        Vec dydt(n_interior);
        const double dx2 = dx * dx;

        // Evaluate boundary prescriptions (Dirichlet values or Neumann fluxes).
        // Check std::function directly (it has an explicit operator bool).
        double u_left_val = left_val ? left_val(t) : 0.0;
        double u_right_val = right_val ? right_val(t) : 0.0;

        for (std::size_t i = 0; i < n_interior; ++i) {
            const double u_i = y[i];

            double u_im1; // u_{i-1}
            double u_ip1; // u_{i+1}

            // left neighbor
            if (i == 0) {
                if (left_bc == pde::BCType::Dirichlet) {
                    u_im1 = u_left_val;
                } else {
                    // Neumann: prescribe du/dx = g(t) at left boundary.
                    // Use ghost-point such that (u1 - u_ghost)/(2*dx) = g
                    // => u_ghost = u1 - 2*dx*g
                    const double g = u_left_val;
                    const double u1 = y[0];
                    const double u_ghost = u1 - 2.0 * dx * g;
                    u_im1 = u_ghost;
                }
            } else {
                u_im1 = y[i - 1];
            }

            // right neighbor
            if (i == n_interior - 1) {
                if (right_bc == pde::BCType::Dirichlet) {
                    u_ip1 = u_right_val;
                } else {
                    // Neumann at right: prescribe du/dx = g(t)
                    // Use ghost-point such that (u_ghost - u_{n})/(2*dx) = g
                    // => u_ghost = u_n + 2*dx*g
                    const double g = u_right_val;
                    const double u_n = y[n_interior - 1];
                    const double u_ghost = u_n + 2.0 * dx * g;
                    u_ip1 = u_ghost;
                }
            } else {
                u_ip1 = y[i + 1];
            }

            dydt[i] = alpha * (u_im1 - 2.0 * u_i + u_ip1) / dx2;
        }

        return dydt;
    };
}
