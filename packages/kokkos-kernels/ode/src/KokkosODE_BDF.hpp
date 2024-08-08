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

#ifndef KOKKOSODE_BDF_HPP
#define KOKKOSODE_BDF_HPP

/// \author Luc Berger-Vergiat (lberge@sandia.gov)
/// \file KokkosODE_BDF.hpp

#include "Kokkos_Core.hpp"
#include "KokkosODE_Types.hpp"
#include "KokkosODE_RungeKutta.hpp"

#include "KokkosODE_BDF_impl.hpp"

namespace KokkosODE {
namespace Experimental {

enum BDF_type : int { BDF1 = 0, BDF2 = 1, BDF3 = 2, BDF4 = 3, BDF5 = 4, BDF6 = 5 };

template <BDF_type T>
struct BDF_coeff_helper {
  using table_type = void;

  BDF_coeff_helper() = default;
};

template <>
struct BDF_coeff_helper<BDF_type::BDF1> {
  using table_type = KokkosODE::Impl::BDF_table<1>;

  BDF_coeff_helper() = default;
};

template <>
struct BDF_coeff_helper<BDF_type::BDF2> {
  using table_type = KokkosODE::Impl::BDF_table<2>;

  BDF_coeff_helper() = default;
};

template <>
struct BDF_coeff_helper<BDF_type::BDF3> {
  using table_type = KokkosODE::Impl::BDF_table<3>;

  BDF_coeff_helper() = default;
};

template <>
struct BDF_coeff_helper<BDF_type::BDF4> {
  using table_type = KokkosODE::Impl::BDF_table<4>;

  BDF_coeff_helper() = default;
};

template <>
struct BDF_coeff_helper<BDF_type::BDF5> {
  using table_type = KokkosODE::Impl::BDF_table<5>;

  BDF_coeff_helper() = default;
};

template <>
struct BDF_coeff_helper<BDF_type::BDF6> {
  using table_type = KokkosODE::Impl::BDF_table<6>;

  BDF_coeff_helper() = default;
};

template <BDF_type T>
struct BDF {
  using table_type = typename BDF_coeff_helper<T>::table_type;

  template <class ode_type, class vec_type, class mv_type, class mat_type, class scalar_type>
  KOKKOS_FUNCTION static void Solve(const ode_type& ode, const scalar_type t_start, const scalar_type t_end,
                                    const int num_steps, const vec_type& y0, const vec_type& y, const vec_type& rhs,
                                    const vec_type& update, const vec_type& scale, const mv_type& y_vecs,
                                    const mv_type& kstack, const mat_type& temp, const mat_type& jac) {
    const table_type table{};

    const double dt = (t_end - t_start) / num_steps;
    double t        = t_start;

    // Load y0 into y_vecs(:, 0)
    for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
      y_vecs(eqIdx, 0) = y0(eqIdx);
    }

    // Compute initial start-up history vectors
    // Using a non adaptive explicit method.
    const int init_steps = table.order - 1;
    if (num_steps < init_steps) {
      return;
    }
    KokkosODE::Experimental::ODE_params params(table.order - 1);
    for (int stepIdx = 0; stepIdx < init_steps; ++stepIdx) {
      KokkosODE::Experimental::RungeKutta<RK_type::RKF45>::Solve(ode, params, t, t + dt, y0, y, update, kstack);

      for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
        y_vecs(eqIdx, stepIdx + 1) = y(eqIdx);
        y0(eqIdx)                  = y(eqIdx);
      }
      t += dt;
    }

    for (int stepIdx = init_steps; stepIdx < num_steps; ++stepIdx) {
      KokkosODE::Impl::BDFStep(ode, table, t, dt, y0, y, rhs, update, scale, y_vecs, temp, jac);

      // Update history
      for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
        y0(eqIdx) = y(eqIdx);
        for (int orderIdx = 0; orderIdx < table.order - 1; ++orderIdx) {
          y_vecs(eqIdx, orderIdx) = y_vecs(eqIdx, orderIdx + 1);
        }
        y_vecs(eqIdx, table.order - 1) = y(eqIdx);
      }
      t += dt;
    }
  }  // Solve()
};

/// \brief BDF Solve integrates an ordinary differential equation
/// using an order and time adaptive BDF method.
///
/// The integration starts with a BDF1 method and adaptively increases
/// or decreases both dt and the order of integration based on error
/// estimators. This function is marked as KOKKOS_FUNCTION so it can
/// be called on host and device.
///
/// \tparam ode_type the type of the ode object to integrated
/// \tparam mv_type a rank-2 view
/// \tparam vec_type a rank-1 view
///
/// \param ode [in]: the ode to integrate
/// \param t_start [in]: time at which the integration starts
/// \param t_end [in]: time at which the integration stops
/// \param initial_step [in]: initial value for dt
/// \param max_step [in]: maximum value for dt
/// \param y0 [in/out]: vector of initial conditions, set to the solution
/// at the end of the integration
/// \param y_new [out]: vector of solution at t_end
/// \param temp [in]: vectors for temporary storage
/// \param temp2 [in]: vectors for temporary storage
template <class ode_type, class mat_type, class vec_type, class scalar_type>
KOKKOS_FUNCTION void BDFSolve(const ode_type& ode, const scalar_type t_start, const scalar_type t_end,
                              const scalar_type initial_step, const scalar_type max_step, const vec_type& y0,
                              const vec_type& y_new, mat_type& temp, mat_type& temp2) {
  using KAT = Kokkos::ArithTraits<scalar_type>;

  // This needs to go away and be pulled out of temp instead...
  auto rhs    = Kokkos::subview(temp, Kokkos::ALL(), 0);
  auto update = Kokkos::subview(temp, Kokkos::ALL(), 1);
  // vec_type rhs("rhs", ode.neqs), update("update", ode.neqs);
  (void)max_step;

  int order = 1, num_equal_steps = 0;
  constexpr scalar_type min_factor = 0.2;
  scalar_type dt                   = initial_step;
  scalar_type t                    = t_start;

  constexpr int max_newton_iters = 10;
  scalar_type atol = 1.0e-6, rtol = 1.0e-3;

  // Compute rhs = f(t_start, y0)
  ode.evaluate_function(t_start, 0, y0, rhs);

  // Check if we need to compute the initial
  // time step size.
  if (initial_step == KAT::zero()) {
    KokkosODE::Impl::initial_step_size(ode, order, t_start, atol, rtol, y0, rhs, temp, dt);
  }

  // Initialize D(:, 0) = y0 and D(:, 1) = dt*rhs
  auto D = Kokkos::subview(temp, Kokkos::ALL(), Kokkos::pair<int, int>(2, 10));
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    D(eqIdx, 0) = y0(eqIdx);
    D(eqIdx, 1) = dt * rhs(eqIdx);
    rhs(eqIdx)  = 0;
  }

  // Now we loop over the time interval [t_start, t_end]
  // and solve our ODE.
  while (t < t_end) {
    KokkosODE::Impl::BDFStep(ode, t, dt, t_end, order, num_equal_steps, max_newton_iters, atol, rtol, min_factor, y0,
                             y_new, rhs, update, temp, temp2);

    for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
      y0(eqIdx) = y_new(eqIdx);
    }
    // printf("t=%f, dt=%f, y={%f, %f, %f}\n", t, dt, y0(0), y0(1), y0(2));
  }
}  // BDFSolve

}  // namespace Experimental
}  // namespace KokkosODE

#endif  // KOKKOSODE_BDF_HPP
