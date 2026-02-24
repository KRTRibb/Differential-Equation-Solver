#pragma once
#include "types.hpp"
#include <functional>
#include <cstddef>

namespace diffeq {
namespace pde {

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
RHS make_heat_rhs_1d(double alpha,
                  std::size_t n_interior,
                  double dx,
                  BCType left_bc = BCType::Dirichlet,
                  BCType right_bc = BCType::Dirichlet,
                  BCFunc left_val = BCFunc(),
                  BCFunc right_val = BCFunc());

RHS make_heat_rhs_2d(double alpha,
                     std::size_t nx, std::size_t ny,
                     double dx, double dy,
                     BCType left_bc = BCType::Dirichlet,
                     BCType right_bc = BCType::Dirichlet,
                     BCType top_bc = BCType::Dirichlet,
                     BCType bottom_bc = BCType::Dirichlet,
                     BCFunc left_val = BCFunc(),
                     BCFunc right_val = BCFunc(),
                     BCFunc top_val = BCFunc(),
                     BCFunc bottom_val = BCFunc());

RHS make_heat_rhs_3d(double alpha,
                     std::size_t nx, std::size_t ny, std::size_t nz,
                     double dx, double dy, double dz,
                     BCType left_bc = BCType::Dirichlet,
                     BCType right_bc = BCType::Dirichlet,
                     BCType top_bc = BCType::Dirichlet,
                     BCType bottom_bc = BCType::Dirichlet,
                     BCType front_bc = BCType::Dirichlet,
                     BCType back_bc = BCType::Dirichlet,
                     BCFunc left_val = BCFunc(),
                     BCFunc right_val = BCFunc(),
                     BCFunc top_val = BCFunc(),
                     BCFunc bottom_val = BCFunc(),
                     BCFunc front_val = BCFunc(),
                     BCFunc back_val = BCFunc());

} // namespace pde
} // namespace diffeq
