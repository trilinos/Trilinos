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

#ifndef KOKKOSODE_NEWTON_IMPL_HPP
#define KOKKOSODE_NEWTON_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_Gesv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"

#include "KokkosODE_Types.hpp"

namespace KokkosODE {
namespace Impl {

template <class system_type, class mat_type, class ini_vec_type, class rhs_vec_type, class update_type,
          class scale_type>
KOKKOS_FUNCTION KokkosODE::Experimental::newton_solver_status NewtonSolve(
    system_type& sys, const KokkosODE::Experimental::Newton_params& params, mat_type& J, mat_type& tmp,
    ini_vec_type& y0, rhs_vec_type& rhs, update_type& update, const scale_type& scale) {
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using value_type           = typename ini_vec_type::non_const_value_type;

  // Define the type returned by nrm2 to store
  // the norm of the residual.
  using norm_type =
      typename Kokkos::Details::InnerProductSpaceTraits<typename ini_vec_type::non_const_value_type>::mag_type;
  sys.residual(y0, rhs);
  const norm_type norm0 = KokkosBlas::serial_nrm2(rhs);
  norm_type norm        = Kokkos::ArithTraits<norm_type>::zero();
  norm_type norm_old    = Kokkos::ArithTraits<norm_type>::zero();
  norm_type norm_new    = Kokkos::ArithTraits<norm_type>::zero();
  norm_type rate        = Kokkos::ArithTraits<norm_type>::zero();

  const norm_type tol = Kokkos::max(10 * Kokkos::ArithTraits<norm_type>::eps() / params.rel_tol,
                                    Kokkos::min(0.03, Kokkos::sqrt(params.rel_tol)));

  // LBV - 07/24/2023: for now assume that we take
  // a full Newton step. Eventually this value can
  // be computed using a line search algorithm to
  // improve convergence for difficult problems.
  const value_type alpha = Kokkos::ArithTraits<value_type>::one();

  // Iterate until maxIts or the tolerance is reached
  for (int it = 0; it < params.max_iters; ++it) {  // handle.maxIters; ++it) {
    // compute initial rhs
    sys.residual(y0, rhs);

    // Solve the following linearized
    // problem at each iteration: J*update=-rhs
    // with J=du/dx, rhs=f(u_n+update)-f(u_n)

    // compute LHS
    sys.jacobian(y0, J);

    // solve linear problem
    int linSolverStat = KokkosBatched::SerialGesv<KokkosBatched::Gesv::StaticPivoting>::invoke(J, update, rhs, tmp);
    KokkosBlas::SerialScale::invoke(-1, update);

    // update solution // x = x + alpha*update
    KokkosBlas::serial_axpy(alpha, update, y0);
    norm = KokkosBlas::serial_nrm2(rhs);

    // Compute rms norm of the scaled update
    for (int idx = 0; idx < sys.neqs; ++idx) {
      norm_new = (update(idx) * update(idx)) / (scale(idx) * scale(idx));
    }
    norm_new = Kokkos::sqrt(norm_new / sys.neqs);
    if ((it > 0) && norm_old > Kokkos::ArithTraits<norm_type>::zero()) {
      rate = norm_new / norm_old;
      if ((rate >= 1) || Kokkos::pow(rate, params.max_iters - it) / (1 - rate) * norm_new > tol) {
        return newton_solver_status::NLS_DIVERGENCE;
      } else if ((norm_new == 0) || ((rate / (1 - rate)) * norm_new < tol)) {
        return newton_solver_status::NLS_SUCCESS;
      }
    }

    if (linSolverStat == 1) {
      Kokkos::printf("NewtonFunctor: Linear solve gesv returned failure! \n");
      return newton_solver_status::LIN_SOLVE_FAIL;
    }

    if ((norm < (params.rel_tol * norm0)) || (it > 0 ? KokkosBlas::serial_nrm2(update) < params.abs_tol : false)) {
      return newton_solver_status::NLS_SUCCESS;
    }

    norm_old = norm_new;
  }
  return newton_solver_status::MAX_ITER;
}

}  // namespace Impl
}  // namespace KokkosODE

#endif  // KOKKOSODE_NEWTON_IMPL_HPP
