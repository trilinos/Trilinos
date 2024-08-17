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

#ifndef KOKKOSBLAS_RUNGEKUTTA_IMPL_HPP
#define KOKKOSBLAS_RUNGEKUTTA_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosODE_RungeKuttaTables_impl.hpp"
#include "KokkosODE_Types.hpp"

namespace KokkosODE {
namespace Impl {

// y_new = y_old + dt*sum(b_i*k_i)    i in [1, nstages]
// k_i = f(t+c_i*dt, y_old+sum(a_{ij}*k_i))  j in [1, i-1]
// we need to compute the k_i and store them as we go
// to use them for k_{i+1} computation.
template <class ode_type, class table_type, class vec_type, class mv_type, class scalar_type>
KOKKOS_FUNCTION void RKStep(ode_type& ode, const table_type& table, const bool adaptivity, scalar_type t,
                            scalar_type dt, const vec_type& y_old, const vec_type& y_new, const vec_type& temp,
                            const mv_type& k_vecs) {
  const int neqs    = ode.neqs;
  const int nstages = table.nstages;

  // first set y_new = y_old
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    y_new(eqIdx) = y_old(eqIdx);
  }

  // now accumulate y_new += dt*b_i*k_i
  {
    // we always start with y_new += dt*b_0*k0
    auto k0 = Kokkos::subview(k_vecs, 0, Kokkos::ALL);
    ode.evaluate_function(t + table.c[0] * dt, dt, y_old, k0);
    for (int eqIdx = 0; eqIdx < neqs; ++eqIdx) {
      y_new(eqIdx) += dt * table.b[0] * k0(eqIdx);
    }
  }

  // Now that we have k0, we can compute all other k_i
  // and accumulate them in y_new.
  for (int stageIdx = 1; stageIdx < nstages; ++stageIdx) {
    for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
      temp(eqIdx) = 0;
    }

    for (int idx = 0; idx < stageIdx; ++idx) {
      for (int eqIdx = 0; eqIdx < neqs; ++eqIdx) {
        temp(eqIdx) += table.a[stageIdx * (stageIdx + 1) / 2 + idx] * k_vecs(idx, eqIdx);
      }
    }
    KokkosBlas::SerialScale::invoke(dt, temp);
    KokkosBlas::serial_axpy(1, y_old, temp);
    auto k = Kokkos::subview(k_vecs, stageIdx, Kokkos::ALL);
    ode.evaluate_function(t + table.c[stageIdx] * dt, dt, temp, k);
    for (int eqIdx = 0; eqIdx < neqs; ++eqIdx) {
      y_new(eqIdx) += dt * table.b[stageIdx] * k(eqIdx);
    }
  }

  // Compute estimation of the error using k_vecs and table.e
  if (adaptivity == true) {
    for (int eqIdx = 0; eqIdx < neqs; ++eqIdx) {
      temp(eqIdx) = 0;
      for (int stageIdx = 0; stageIdx < nstages; ++stageIdx) {
        temp(eqIdx) += dt * table.e[stageIdx] * k_vecs(stageIdx, eqIdx);
      }
    }
  }
}  // RKStep

template <class ode_type, class table_type, class vec_type, class mv_type, class scalar_type>
KOKKOS_FUNCTION Experimental::ode_solver_status RKSolve(const ode_type& ode, const table_type& table,
                                                        const KokkosODE::Experimental::ODE_params& params,
                                                        const scalar_type t_start, const scalar_type t_end,
                                                        const vec_type& y0, const vec_type& y, const vec_type& temp,
                                                        const mv_type& k_vecs) {
  constexpr scalar_type error_threshold = 1;
  bool adapt                            = params.adaptivity;
  bool dt_was_reduced;
  if (std::is_same_v<table_type, ButcherTableau<0, 0>>) {
    adapt = false;
  }

  // Set current time and initial time step
  scalar_type t_now = t_start;
  scalar_type dt    = (t_end - t_start) / params.max_steps;

  // Loop over time steps to integrate ODE
  for (int stepIdx = 0; (stepIdx < params.max_steps) && (t_now <= t_end); ++stepIdx) {
    // Check that the step attempted is not putting
    // the solution past t_end, otherwise shrink dt
    if (t_end < t_now + dt) {
      dt = t_end - t_now;
    }

    // Set error to be arbitrarily larger than our threshold
    // so we can pass the initial check. Also reset
    // dt_was_reduced to false for current time step.
    scalar_type error = 2 * error_threshold;
    scalar_type tol   = 0;
    dt_was_reduced    = false;

    // Take tentative steps until the requested error
    // is met. This of course only works for adaptive
    // solvers, for fix time steps we simply do not
    // compute and check what error of the current step
    while (error_threshold < error) {
      // Take a step of Runge-Kutta integrator
      RKStep(ode, table, adapt, t_now, dt, y0, y, temp, k_vecs);

      // Compute the largest error and decide on
      // the size of the next time step to take.
      error = 0;
      if (adapt) {
        // Compute the error
        for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
          error = Kokkos::max(error, Kokkos::abs(temp(eqIdx)));
          tol   = Kokkos::max(
              tol, params.abs_tol + params.rel_tol * Kokkos::max(Kokkos::abs(y(eqIdx)), Kokkos::abs(y0(eqIdx))));
        }
        error = error / tol;

        // Reduce the time step if error
        // is too large and current step
        // is rejected.
        if (error > 1) {
          dt             = dt * Kokkos::max(0.2, 0.8 / Kokkos::pow(error, 1 / table.order));
          dt_was_reduced = true;
        }

        if (dt < params.min_step_size) return Experimental::ode_solver_status::MIN_SIZE;
      }
    }

    // Update time and initial condition for next time step
    t_now += dt;
    for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
      y0(eqIdx) = y(eqIdx);
    }

    if (t_now < t_end) {
      if (adapt && !dt_was_reduced && error < 0.5) {
        // Compute new time increment
        dt = dt * Kokkos::min(10.0, Kokkos::max(2.0, 0.9 * Kokkos::pow(error, 1 / table.order)));
      }
    } else {
      return Experimental::ode_solver_status::SUCCESS;
    }
  }

  if (t_now < t_end) return Experimental::ode_solver_status::MAX_STEP;

  return Experimental::ode_solver_status::SUCCESS;
}  // RKSolve

}  // namespace Impl
}  // namespace KokkosODE

#endif  // KOKKOSBLAS_RUNGEKUTTA_IMPL_HPP
