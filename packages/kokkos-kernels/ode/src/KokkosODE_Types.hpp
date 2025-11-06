// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSODE_TYPES_HPP
#define KOKKOSODE_TYPES_HPP

namespace KokkosODE {
namespace Experimental {

enum ode_solver_status { SUCCESS = 0, MAX_STEP = 1, MIN_SIZE = 2 };

struct ODE_params {
  bool adaptivity;
  int num_steps, max_steps;
  double abs_tol, rel_tol, min_step_size;

  KOKKOS_FUNCTION
  ODE_params()
      : adaptivity(true), num_steps(100), max_steps(10000), abs_tol(1e-12), rel_tol(1e-6), min_step_size(1e-9) {}

  // Constructor that only specify the desired number of steps.
  // In this case no adaptivity is provided, the time step will
  // be constant such that dt = (tend - tstart) / num_steps;
  KOKKOS_FUNCTION
  ODE_params(const int num_steps_)
      : adaptivity(false),
        num_steps(num_steps_),
        max_steps(num_steps_ + 1),
        abs_tol(1e-12),
        rel_tol(1e-6),
        min_step_size(0) {}

  /// ODE_parms construtor for adaptive time stepping.
  KOKKOS_FUNCTION
  ODE_params(const int num_steps_, const int max_steps_, const double abs_tol_, const double rel_tol_,
             const double min_step_size_)
      : adaptivity(true),
        num_steps(num_steps_),
        max_steps(max_steps_),
        abs_tol(abs_tol_),
        rel_tol(rel_tol_),
        min_step_size(min_step_size_) {}
};

enum newton_solver_status : int {
  NLS_SUCCESS    = 0,
  MAX_ITER       = 1,
  LIN_SOLVE_FAIL = 2,
  NLS_DIVERGENCE = 3,
};

struct Newton_params {
  int max_iters, iters = 0;
  double abs_tol, rel_tol;

  // Constructor that sets basic solver parameters
  // used while solving the nonlinear system
  // int max_iters_  [in]: maximum number of iterations allowed
  // double abs_tol_ [in]: absolute tolerance to reach for successful solve
  // double rel_tol_ [in]: relative tolerance to reach for successful solve
  KOKKOS_FUNCTION
  Newton_params(const int max_iters_, const double abs_tol_, const double rel_tol_)
      : max_iters(max_iters_), abs_tol(abs_tol_), rel_tol(rel_tol_) {}
};

}  // namespace Experimental
}  // namespace KokkosODE
#endif  // KOKKOSODE_TYPES_HPP
