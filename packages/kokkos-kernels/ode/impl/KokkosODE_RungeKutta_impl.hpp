// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS_RUNGEKUTTA_IMPL_HPP
#define KOKKOSBLAS_RUNGEKUTTA_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosODE_RungeKuttaTables_impl.hpp"
#include "KokkosODE_Types.hpp"

#include "iostream"

namespace KokkosODE {
namespace Impl {

// This algorithm is mostly derived from
// E. Hairer, S. P. Norsett G. Wanner,
// "Solving Ordinary Differential Equations I:
// Nonstiff Problems", Sec. II.4.
// Note that all floating point values below
// have been heuristically selected for
// convergence performance.
template <class ode_type, class mat_type, class vec_type, class res_type, class scalar_type>
KOKKOS_FUNCTION void first_step_size(const ode_type ode, const int order, const scalar_type t0, const scalar_type atol,
                                     const scalar_type rtol, const vec_type& y0, const res_type& f0, const vec_type y1,
                                     const mat_type temp, scalar_type& dt_ini) {
  using KAT = KokkosKernels::ArithTraits<scalar_type>;

  // Extract subviews to store intermediate data
  auto f1 = Kokkos::subview(temp, 1, Kokkos::ALL());

  // Compute norms for y0 and f0
  double n0 = KAT::zero(), n1 = KAT::zero(), dt0, scale;
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    scale = atol + rtol * Kokkos::abs(y0(eqIdx));
    n0 += Kokkos::pow(y0(eqIdx) / scale, 2);
    n1 += Kokkos::pow(f0(eqIdx) / scale, 2);
  }
  n0 = Kokkos::sqrt(n0) / Kokkos::sqrt(ode.neqs);
  n1 = Kokkos::sqrt(n1) / Kokkos::sqrt(ode.neqs);

  // Select dt0
  if ((n0 < 1e-5) || (n1 < 1e-5)) {
    dt0 = 1e-6;
  } else {
    dt0 = 0.01 * n0 / n1;
  }

  // Estimate y at t0 + dt0
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    y1(eqIdx) = y0(eqIdx) + dt0 * f0(eqIdx);
  }

  // Compute f at t0+dt0 and y1,
  // then compute the norm of f(t0+dt0, y1) - f(t0, y0)
  scalar_type n2 = KAT::zero();
  ode.evaluate_function(t0 + dt0, dt0, y1, f1);
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    n2 += Kokkos::pow((f1(eqIdx) - f0(eqIdx)) / (atol + rtol * Kokkos::abs(y0(eqIdx))), 2);
  }
  n2 = Kokkos::sqrt(n2) / (dt0 * Kokkos::sqrt(ode.neqs));

  // Finally select initial time step dt_ini
  if ((n1 <= 1e-15) && (n2 <= 1e-15)) {
    dt_ini = Kokkos::max(1e-6, dt0 * 1e-3);
  } else {
    dt_ini = Kokkos::pow(0.01 / Kokkos::max(n1, n2), KAT::one() / order);
  }

  dt_ini = Kokkos::min(100 * dt0, dt_ini);

  // Zero out temp variables just to be safe...
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    f0(eqIdx) = 0.0;
    y1(eqIdx) = 0.0;
    f1(eqIdx) = 0.0;
  }
}  // first_step_size

// y_new = y_old + dt*sum(b_i*k_i)    i in [1, nstages]
// k_i = f(t+c_i*dt, y_old+sum(a_{ij}*k_i))  j in [1, i-1]
// we need to compute the k_i and store them as we go
// to use them for k_{i+1} computation.
template <class ode_type, class table_type, class vec_type, class mv_type, class scalar_type>
KOKKOS_FUNCTION void RKStep(ode_type& ode, const table_type& table, scalar_type t, scalar_type dt,
                            const vec_type& y_old, const vec_type& y_new, const vec_type& temp, const mv_type& k_vecs) {
  const int neqs        = ode.neqs;
  constexpr int nstages = table_type::nstages;

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
}  // RKStep

// Note that the control values for
// time step increase/decrease are
// heuristically chosen based on
// L. F. Shampine and M. W. Reichelt
// "The Matlab ODE suite" SIAM J. Sci.
// Comput. Vol. 18, No. 1, pp. 1-22
// Jan. 1997
template <class ode_type, class table_type, class vec_type, class mv_type, class scalar_type>
KOKKOS_FUNCTION Experimental::ode_solver_status RKSolve(const ode_type& ode, const table_type& table,
                                                        const KokkosODE::Experimental::ODE_params& params,
                                                        const scalar_type t_start, const scalar_type t_end,
                                                        const vec_type& y0, const vec_type& y, const vec_type& temp,
                                                        const mv_type& k_vecs, int* const step_count) {
  constexpr scalar_type error_threshold = 1;
  scalar_type error_n;
  bool adapt = params.adaptivity;
  bool dt_was_reduced;
  if constexpr (std::is_same_v<table_type, ButcherTableau<0, 0>>) {
    adapt = false;
  }

  // Set current time and initial time step
  scalar_type t_now = t_start, dt = 0.0;
  if (adapt == true) {
    ode.evaluate_function(t_start, 0, y0, temp);
    first_step_size(ode, table_type::order, t_start, params.abs_tol, params.rel_tol, y0, temp, y, k_vecs, dt);
    if (dt < params.min_step_size) {
      dt = params.min_step_size;
    }
  } else {
    dt = (t_end - t_start) / params.num_steps;
  }

  *step_count = 0;

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
    // compute and check the error of the current step
    while (error_threshold < error) {
      // Take a step of Runge-Kutta integrator
      RKStep(ode, table, t_now, dt, y0, y, temp, k_vecs);

      // Compute the largest error and decide on
      // the size of the next time step to take.
      error = 0;

      // Compute estimation of the error using k_vecs and table.e
      if (adapt == true) {
        // Compute the error
        for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
          tol     = params.abs_tol + params.rel_tol * Kokkos::max(Kokkos::abs(y(eqIdx)), Kokkos::abs(y0(eqIdx)));
          error_n = 0;
          for (int stageIdx = 0; stageIdx < table.nstages; ++stageIdx) {
            error_n += dt * table.e[stageIdx] * k_vecs(stageIdx, eqIdx);
          }
          error += (error_n * error_n) / (tol * tol);
        }
        error = Kokkos::sqrt(error / ode.neqs);

        // Reduce the time step if error
        // is too large and current step
        // is rejected.
        if (error > 1) {
          dt             = dt * Kokkos::max(0.2, 0.8 * Kokkos::pow(error, -1.0 / table.order));
          dt_was_reduced = true;
        }

        if (dt < params.min_step_size) return Experimental::ode_solver_status::MIN_SIZE;
      }
    }

    // Update time and initial condition for next time step
    t_now += dt;
    *step_count += 1;
    for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
      y0(eqIdx) = y(eqIdx);
    }

    if (t_now < t_end) {
      if (adapt && !dt_was_reduced && error < 0.5) {
        // Compute new time increment
        dt = dt * Kokkos::min(10.0, Kokkos::max(2.0, 0.9 * Kokkos::pow(error, -1.0 / table.order)));
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
