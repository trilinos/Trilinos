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

#ifndef KOKKOSODE_RUNGEKUTTA_HPP
#define KOKKOSODE_RUNGEKUTTA_HPP

/// \author Luc Berger-Vergiat (lberge@sandia.gov)
/// \file KokkosODE_RungeKutta.hpp

#include "Kokkos_Core.hpp"
#include "KokkosODE_Types.hpp"

#include "KokkosODE_RungeKutta_impl.hpp"

namespace KokkosODE {
namespace Experimental {

/// \brief RK_type is an enum tye that conveniently
/// describes the Runge-Kutta methods implemented.
enum RK_type : int {
  RKFE  = 0,  ///< Forward Euler method (no adaptivity available for this method)
  RKEH  = 1,  ///< Euler-Heun method
  RKF12 = 2,  ///< Fehlberg order 2 method
  RKBS  = 3,  ///< Bogacki-Shampine method
  RK4   = 4,  ///< Runge-Kutta classic order 4 method
  RKF45 = 5,  ///< Fehlberg order 5 method
  RKCK  = 6,  ///< Cash-Karp method
  RKDP  = 7   ///< Dormand-Prince method
};

template <RK_type T>
struct RK_Tableau_helper {
  using table_type = void;
};

template <>
struct RK_Tableau_helper<RK_type::RKFE> {
  using table_type = KokkosODE::Impl::ButcherTableau<0, 0>;
};

template <>
struct RK_Tableau_helper<RK_type::RKEH> {
  using table_type = KokkosODE::Impl::ButcherTableau<1, 1>;
};

template <>
struct RK_Tableau_helper<RK_type::RKF12> {
  using table_type = KokkosODE::Impl::ButcherTableau<1, 2>;
};

template <>
struct RK_Tableau_helper<RK_type::RKBS> {
  using table_type = KokkosODE::Impl::ButcherTableau<2, 3>;
};

template <>
struct RK_Tableau_helper<RK_type::RK4> {
  using table_type = KokkosODE::Impl::ButcherTableau<3, 3>;
};

template <>
struct RK_Tableau_helper<RK_type::RKF45> {
  using table_type = KokkosODE::Impl::ButcherTableau<4, 5>;
};

template <>
struct RK_Tableau_helper<RK_type::RKCK> {
  using table_type = KokkosODE::Impl::ButcherTableau<4, 5, 1>;
};

template <>
struct RK_Tableau_helper<RK_type::RKDP> {
  using table_type = KokkosODE::Impl::ButcherTableau<4, 6>;
};

/// \brief Unspecialized version of the RungeKutta solvers
///
/// \tparam RK_type an RK_type enum value used to specify
///         which Runge Kutta method is to be used.
template <RK_type T>
struct RungeKutta {
  using table_type = typename RK_Tableau_helper<T>::table_type;

  /// \brief order returns the convergence order of the method
  KOKKOS_FUNCTION
  static int order() { return table_type::order; }

  /// \brief num_stages returns the number of stages used by the method
  KOKKOS_FUNCTION
  static int num_stages() { return table_type::nstages; }

  /// \brief Solve integrates an ordinary differential equation
  ///
  /// The integration is carried with the method specified as template
  /// parameter to the RungeKutta struct. This method is static and
  /// marked as KOKKOS_FUNCTION so it can be used on host and device.
  ///
  /// \tparam ode_type the type of the ode object to integrated
  /// \tparam vec_type a rank-1 view
  /// \tparam mv_type a rank-2 view
  /// \tparam scalar_type a floating point type
  ///
  /// \param ode [in]: the ode to integrate
  /// \param params [in]: standard input parameters of ODE integrators
  /// \param t_start [in]: time at which the integration starts
  /// \param t_end [in]: time at which the integration stops
  /// \param y0 [in/out]: vector of initial conditions, set to the solution
  /// at the end of the integration
  /// \param y [out]: vector of solution at t_end
  /// \param temp [in]: vector for temporary storage
  /// \param k_vecs [in]: vectors for temporary storage
  ///
  /// \return ode_solver_status an enum that describes success of failure
  /// of the integration method once it at terminated.
  template <class ode_type, class vec_type, class mv_type, class scalar_type>
  KOKKOS_FUNCTION static ode_solver_status Solve(const ode_type& ode, const KokkosODE::Experimental::ODE_params& params,
                                                 const scalar_type t_start, const scalar_type t_end, const vec_type& y0,
                                                 const vec_type& y, const vec_type& temp, const mv_type& k_vecs) {
    table_type table;
    return KokkosODE::Impl::RKSolve(ode, table, params, t_start, t_end, y0, y, temp, k_vecs);
  }
};

}  // namespace Experimental
}  // namespace KokkosODE
#endif  // KOKKOSODE_RUNGEKUTTA_HPP
