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
template <class ode_type, class table_type, class vec_type, class mv_type,
          class scalar_type>
KOKKOS_FUNCTION void RKStep(ode_type& ode, const table_type& table,
                            const bool adaptivity, scalar_type t,
                            scalar_type dt, const vec_type& y_old,
                            const vec_type& y_new, const vec_type& temp,
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
    auto k0 = Kokkos::subview(k_vecs, Kokkos::ALL, 0);
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
        temp(eqIdx) +=
            table.a[stageIdx * (stageIdx + 1) / 2 + idx] * k_vecs(eqIdx, idx);
      }
    }
    KokkosBlas::SerialScale::invoke(dt, temp);
    KokkosBlas::serial_axpy(1, y_old, temp);
    auto k = Kokkos::subview(k_vecs, Kokkos::ALL, stageIdx);
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
        temp(eqIdx) += dt * table.e[stageIdx] * k_vecs(eqIdx, stageIdx);
      }
    }
  }
}  // RKStep

template <class ode_type, class table_type, class vec_type, class mv_type,
          class scalar_type>
KOKKOS_FUNCTION Experimental::ode_solver_status RKSolve(
    const ode_type& ode, const table_type& table,
    const KokkosODE::Experimental::ODE_params& params,
    const scalar_type t_start, const scalar_type t_end, const vec_type& y0,
    const vec_type& y, const vec_type& temp, const mv_type& k_vecs) {
  constexpr scalar_type error_threshold = 1;
  bool adapt                            = params.adaptivity;
  if (std::is_same_v<table_type, ButcherTableau<0, 0>>) {
    adapt = false;
  }

  scalar_type dt = (t_end - t_start) / params.max_steps;
  scalar_type t  = t_start;
  for (int stepIdx = 0; (stepIdx < params.max_steps) && (t < t_end);
       ++stepIdx) {
    // Set err to be arbitrarily larger than our threshold of 1
    scalar_type error = 2 * error_threshold;
    scalar_type tol   = 0;
    while (error_threshold < error) {
      // Take a step of Runge-Kutta integrator
      RKStep(ode, table, adapt, t, dt, y0, y, temp, k_vecs);

      // Compute the largest error and decide on
      // the size of the next time step to take.
      error = 0;
      if (adapt) {
        // Compute the error
        for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
          error = Kokkos::max(error, Kokkos::abs(temp(eqIdx)));
          tol   = Kokkos::max(
              tol, params.abs_tol +
                       params.rel_tol * Kokkos::max(Kokkos::abs(y(eqIdx)),
                                                    Kokkos::abs(y0(eqIdx))));
        }
        error = error / tol;

        // Reduce the time step if error
        // is too large and current step
        // is rejected.
        if (error > 1) {
          dt = dt * Kokkos::max(0.2, 0.8 / Kokkos::pow(error, 1 / table.order));
        }
        if (dt < params.min_step_size)
          return Experimental::ode_solver_status::MIN_SIZE;
      }
    }

    // Update y0 to stage the next time step.
    for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
      y0(eqIdx) = y(eqIdx);
    }

    if (t < t_end) {
      // We may want to print the evolution of the solution over time
      // with something similar to the statement below but will need
      // to generalize it and make it GPU friendly first, also it
      // should be guarded when not doing a debug run, this prints
      // a lot...
      // std::cout << " step " << stepIdx << " t=" << t << ", y={";
      // for(int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
      //   std::cout << y(eqIdx) << " ";
      // }
      // std::cout << "}" << std::endl;
      if (adapt) {
        // Compute new time increment
        dt = dt *
             Kokkos::min(
                 10.0,
                 Kokkos::max(2.0, 0.9 * Kokkos::pow(error, 1 / table.order)));
      } else {
        // Use same increment
        t += dt;
      }
    } else {
      return Experimental::ode_solver_status::SUCCESS;
    }
  }

  if (t < t_end) return Experimental::ode_solver_status::MAX_STEP;

  return Experimental::ode_solver_status::SUCCESS;
}  // RKSolve

}  // namespace Impl
}  // namespace KokkosODE

#endif  // KOKKOSBLAS_RUNGEKUTTA_IMPL_HPP
