// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSODE_NEWTON_IMPL_HPP
#define KOKKOSODE_NEWTON_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_Gesv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosKernels_InnerProductSpaceTraits.hpp"

#include "KokkosBatched_Getrf.hpp"
#include "KokkosBatched_Getrs.hpp"
#include "KokkosODE_Types.hpp"

namespace KokkosODE {
namespace Impl {

template <class system_type, class mat_type, class ini_vec_type, class rhs_vec_type, class update_type,
          class scale_type>
KOKKOS_FUNCTION KokkosODE::Experimental::newton_solver_status NewtonSolve(
    system_type& sys, const KokkosODE::Experimental::Newton_params& params, mat_type& J, mat_type& tmp,
    ini_vec_type& y0, rhs_vec_type& rhs, update_type& update, const scale_type& scale, int& newton_iterations) {
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using value_type           = typename ini_vec_type::non_const_value_type;

  // Define the type returned by nrm2 to store
  // the norm of the residual.
  using norm_type =
      typename KokkosKernels::Details::InnerProductSpaceTraits<typename ini_vec_type::non_const_value_type>::mag_type;
  sys.residual(y0, rhs);
  const norm_type norm0 = KokkosBlas::serial_nrm2(rhs);
  norm_type norm        = KokkosKernels::ArithTraits<norm_type>::zero();
  norm_type norm_old    = KokkosKernels::ArithTraits<norm_type>::zero();
  norm_type rate        = KokkosKernels::ArithTraits<norm_type>::zero();

  const norm_type tol = Kokkos::max(10 * KokkosKernels::ArithTraits<norm_type>::eps() / params.rel_tol,
                                    Kokkos::min(0.03, Kokkos::sqrt(params.rel_tol)));

  // LBV - 07/24/2023: for now assume that we take
  // a full Newton step. Eventually this value can
  // be computed using a line search algorithm to
  // improve convergence for difficult problems.
  const value_type alpha = KokkosKernels::ArithTraits<value_type>::one();

  // Iterate until maxIts or the tolerance is reached
  newton_iterations = 0;
  for (int it = 0; it < params.max_iters; ++it) {  // handle.maxIters; ++it) {
    newton_iterations++;

    // compute initial rhs
    sys.residual(y0, rhs);

    // Solve the following linearized
    // problem at each iteration: J*update=-rhs
    // with J=du/dx, rhs=f(u_n+update)-f(u_n)

    // compute LHS
    sys.jacobian(y0, J);

    // Solve the linear problem J*update = rhs using LU with partial
    // (row) pivoting. The static pivoting implemented in
    // KokkosBatched::SerialGesv uses a greedy row-to-column assignment
    // that can fail spuriously on well-conditioned sparse systems
    // (e.g. chemistry Jacobians with one dense row), so it is only kept
    // as a fallback for systems larger than the stack pivot buffer.
    int linSolverStat    = 0;
    using mat_value_type = typename mat_type::non_const_value_type;
    if (sys.neqs * sizeof(int) <= tmp.size() * sizeof(mat_value_type)) {
      Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> piv(
          reinterpret_cast<int*>(tmp.data()), sys.neqs);
      for (int idx = 0; idx < sys.neqs; ++idx) update(idx) = rhs(idx);
      linSolverStat = KokkosBatched::SerialGetrf<KokkosBatched::Algo::Getrf::Unblocked>::invoke(J, piv);
      if (linSolverStat == 0) {
        linSolverStat = KokkosBatched::SerialGetrs<KokkosBatched::Trans::NoTranspose,
                                                   KokkosBatched::Algo::Getrs::Unblocked>::invoke(J, piv, update);
      }
    } else {
      linSolverStat = KokkosBatched::SerialGesv<KokkosBatched::Gesv::StaticPivoting>::invoke(J, update, rhs, tmp);
    }

    // Return before touching y0 if the linear solve failed: applying the
    // (stale or partial) update would corrupt the iterate.
    if (linSolverStat != 0) {
      return newton_solver_status::LIN_SOLVE_FAIL;
    }

    KokkosBlas::SerialScale::invoke(-1, update);

    // Compute the rms norm of the scaled update and check for divergence
    // before applying the update to y0
    norm_type norm_new = KokkosKernels::ArithTraits<norm_type>::zero();
    for (int idx = 0; idx < sys.neqs; ++idx) {
      norm_new += (update(idx) * update(idx)) / (scale(idx) * scale(idx));
    }
    norm_new = Kokkos::sqrt(norm_new / sys.neqs);
    if ((it > 0) && norm_old > KokkosKernels::ArithTraits<norm_type>::zero()) {
      rate = norm_new / norm_old;
      if ((rate >= 1) || Kokkos::pow(rate, params.max_iters - it) / (1 - rate) * norm_new > tol) {
        return newton_solver_status::NLS_DIVERGENCE;
      }
    }

    // update solution // x = x + alpha*update
    KokkosBlas::serial_axpy(alpha, update, y0);
    norm = KokkosBlas::serial_nrm2(rhs);

    if ((it > 0) && norm_old > KokkosKernels::ArithTraits<norm_type>::zero() &&
        ((norm_new == 0) || ((rate / (1 - rate)) * norm_new < tol))) {
      return newton_solver_status::NLS_SUCCESS;
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
