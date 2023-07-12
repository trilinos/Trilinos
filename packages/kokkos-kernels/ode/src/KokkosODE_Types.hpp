//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSODE_TYPES_HPP
#define KOKKOSODE_TYPES_HPP

namespace KokkosODE {
namespace Experimental {

enum ode_solver_status { SUCCESS = 0, MAX_STEP = 1, MIN_SIZE = 2 };

struct ODE_params {
  bool adaptivity;
  int num_steps, max_steps;
  double abs_tol, rel_tol, min_step_size;

  // Constructor that only specify the desired number of steps.
  // In this case no adaptivity is provided, the time step will
  // be constant such that dt = (tend - tstart) / num_steps;
  KOKKOS_FUNCTION
  ODE_params(const int num_steps_)
      : adaptivity(false),
        num_steps(num_steps_),
        max_steps(num_steps_),
        abs_tol(0),
        rel_tol(0),
        min_step_size(0) {}

  /// ODE_parms construtor for adaptive time stepping.
  KOKKOS_FUNCTION
  ODE_params(const int num_steps_, const int max_steps_, const double abs_tol_,
             const double rel_tol_, const double min_step_size_)
      : adaptivity(true),
        num_steps(num_steps_),
        max_steps(max_steps_),
        abs_tol(abs_tol_),
        rel_tol(rel_tol_),
        min_step_size(min_step_size_) {}
};

}  // namespace Experimental
}  // namespace KokkosODE
#endif  // KOKKOSODE_TYPES_HPP
