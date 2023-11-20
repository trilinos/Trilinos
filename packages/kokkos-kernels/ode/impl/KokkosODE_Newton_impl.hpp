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

template <class system_type, class mat_type, class vec_type>
KOKKOS_FUNCTION KokkosODE::Experimental::newton_solver_status NewtonSolve(
    system_type& sys, const KokkosODE::Experimental::Newton_params& params,
    mat_type& J, mat_type& tmp, vec_type& y0, vec_type& rhs, vec_type& update) {
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using value_type           = typename vec_type::non_const_value_type;

  // Define the type returned by nrm2 to store
  // the norm of the residual.
  using norm_type = typename Kokkos::Details::InnerProductSpaceTraits<
      typename vec_type::non_const_value_type>::mag_type;
  norm_type norm = Kokkos::ArithTraits<norm_type>::zero();

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
    norm = KokkosBlas::serial_nrm2(rhs);

    if ((norm < params.rel_tol) ||
        (it > 0 ? KokkosBlas::serial_nrm2(update) < params.abs_tol : false)) {
      return newton_solver_status::NLS_SUCCESS;
    }

    // compute LHS
    sys.jacobian(y0, J);

    // solve linear problem
    int linSolverStat =
        KokkosBatched::SerialGesv<KokkosBatched::Gesv::StaticPivoting>::invoke(
            J, update, rhs, tmp);
    KokkosBlas::SerialScale::invoke(-1, update);

    if (linSolverStat == 1) {
#if KOKKOS_VERSION < 40199
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "NewtonFunctor: Linear solve gesv returned failure! \n");
#else
      Kokkos::printf("NewtonFunctor: Linear solve gesv returned failure! \n");
#endif
      return newton_solver_status::LIN_SOLVE_FAIL;
    }

    // update solution // x = x + alpha*update
    KokkosBlas::serial_axpy(alpha, update, y0);
  }
  return newton_solver_status::MAX_ITER;
}

}  // namespace Impl
}  // namespace KokkosODE

#endif  // KOKKOSODE_NEWTON_IMPL_HPP
