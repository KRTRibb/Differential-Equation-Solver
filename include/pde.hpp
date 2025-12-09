#pragma once
#include "types.hpp"
#include <functional>
#include <cstddef>

namespace diffeq {
namespace pde {

// Simple boundary condition types used by the helpers
enum class BCType { Dirichlet, Neumann };

// Boundary value / flux signature: value at a boundary as a function of time
using BCFunc = std::function<double(double)>;

// Build an RHS for the 1D heat equation using the method of lines.
// The returned RHS expects a state vector `y` containing the values at
// the `n_interior` grid points (interior points only; boundary values
// are supplied via BCs). The spatial spacing `dx` is the distance
// between adjacent nodes (interior-to-interior spacing).
//
// u_t = alpha * u_xx
//
// left_bc/right_bc select Dirichlet (prescribed u) or Neumann
// (prescribed du/dx). For Dirichlet, `left_val`/`right_val` should
// provide u(t) at the boundary. For Neumann, they should provide the
// derivative du/dx at the boundary. If a BCFunc is empty it is treated
// as zero.
RHS make_heat_rhs(double alpha,
                  std::size_t n_interior,
                  double dx,
                  BCType left_bc = BCType::Dirichlet,
                  BCType right_bc = BCType::Dirichlet,
                  BCFunc left_val = BCFunc(),
                  BCFunc right_val = BCFunc());

} // namespace pde
} // namespace diffeq
