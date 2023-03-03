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

#ifndef __KOKKOSBATCHED_ODE_NEWTON_HPP__
#define __KOKKOSBATCHED_ODE_NEWTON_HPP__

#include "Kokkos_Core.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_Gesv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"

namespace KokkosBlas {
namespace Impl {

enum class NewtonSolverStatus { Converged = 0, LinearSolveFailure, MaxIters };

std::ostream& operator<<(std::ostream& os, NewtonSolverStatus& status) {
  switch (status) {
    case NewtonSolverStatus::Converged: os << "Newton Solver Converged!"; break;
    case NewtonSolverStatus::LinearSolveFailure:
      os << "Newton: Linear Solver Failure";
      break;
    case NewtonSolverStatus::MaxIters:
      os << "Newton reached maximum iterations without convergence.";
      break;
  }
  return os;
}

/// \brief NewtonHandle
///
/// This handle is used to pass information between the Newton Solver and
/// the calling code.
///
/// \tparam: NormViewType: Type of view used to store the residual convergence
/// history

template <class NormViewType>
struct NewtonHandle {
  using norm_type = typename NormViewType::non_const_value_type;

  NormViewType lastResidual;  // Residual of last successful iteration
  typename NormViewType::HostMirror lastResidualHost;

  // NormViewType  residual_norms;
  // TODO: Making these public for now. Should make private and access
  // via setters and getters?
  int maxIters;           // Maximum number of Newton steps
  norm_type relativeTol;  // Relative convergence tolerance
  bool debug_mode;        // Returns extra verbose output if true.

  NewtonHandle(int _maxIters = 25, double _relativeTol = 1.0e-6,
               bool _debug = false)
      : lastResidual("ending Residual norm", 1),
        lastResidualHost("end res norm host", 1),
        maxIters(_maxIters),
        relativeTol(_relativeTol),
        debug_mode(_debug) {}

  KOKKOS_FUNCTION
  void set_residual(const norm_type val) const { lastResidual(0) = val; }

  KOKKOS_FUNCTION
  norm_type get_residual() const { return lastResidual(0); }

  norm_type get_residual_host() const {
    Kokkos::deep_copy(lastResidualHost, lastResidual);
    return lastResidualHost(0);
  }

};  // NewtonHandle

/// \brief Newton Functor:
/// Solves the nonlinear system F(x) = 0
/// where F is a map from R^n to R^n.
/// \tparam System: Struct that allows the evaluation
///         of the residual and jacobian using the
///         residual() and jacobian() methods.
/// \tparam Matrix: rank-2 view-type
/// \tparam XVector: rank-1 view-type
/// \tparam YVector: rank-1 view-type
/// \param
/// \param X [in]: Input vector X, a rank 1 view
/// \param Y [in/out]: Output vector Y, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///
template <class System, class Matrix, class XVector, class YVector,
          class NewtonHandleType>
struct NewtonFunctor {
  using execution_space = typename YVector::execution_space;
  using yvalue_type     = typename YVector::non_const_value_type;
  using norm_type       = typename NewtonHandleType::norm_type;

  System sys;
  XVector x;
  YVector rhs;
  NewtonHandleType handle;

  Matrix J, tmp;
  XVector update;

  NewtonFunctor(System _sys, XVector _x, YVector _rhs,
                NewtonHandleType& _handle)
      : sys(_sys), x(_x), rhs(_rhs), handle(_handle) {
    J      = Matrix("Jacobian", x.extent(0), x.extent(0));
    tmp    = Matrix("Jacobian", x.extent(0), x.extent(0) + 4);
    update = XVector("update", x.extent(0));
  }

  KOKKOS_INLINE_FUNCTION
  NewtonSolverStatus solve() const {
    norm_type norm    = Kokkos::ArithTraits<norm_type>::zero();
    yvalue_type alpha = Kokkos::ArithTraits<yvalue_type>::one();
    handle.set_residual(-1);  // init to dummy value

    // Iterate until maxIts or the tolerance is reached
    for (int it = 0; it < handle.maxIters; ++it) {
      // compute initial rhs
      sys.residual(x, rhs);
      if (handle.debug_mode) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("NewtonFunctor: r=");
        for (int k = 0; k < rhs.extent_int(0); k++) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f \n", rhs(k));
        }
      }

      // Solve the following linearized
      // problem at each step: J*update=-rhs
      // with J=du/dx, rhs=f(u_n+update)-f(u_n)
      norm = KokkosBlas::serial_nrm2(rhs);
      handle.set_residual(norm);

      if (handle.debug_mode) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "NewtonFunctor: Iteration: %d  Current res norm is: %e \n Current "
            "soln is:\n",
            it, (double)handle.get_residual());
        for (int k = 0; k < x.extent_int(0); k++) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f \n", x(k));
        }
      }

      if (norm < handle.relativeTol) {
        // Problem solved, exit the functor
        if (handle.debug_mode) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "NewtonFunctor: Newton solver converged! Ending norm is: %e \n "
              "Solution x is: "
              "\n",
              norm);
          for (int k = 0; k < x.extent_int(0); k++) {
            KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f \n", x(k));
          }
        }
        return NewtonSolverStatus::Converged;
      }

      // compute LHS
      sys.jacobian(x, J);

      // solve linear problem
      int linSolverStat = KokkosBatched::SerialGesv<
          KokkosBatched::Gesv::StaticPivoting>::invoke(J, update, rhs, tmp);
      KokkosBlas::SerialScale::invoke(-1, update);

      if (handle.debug_mode) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "NewtonFunctor: Print linear solve solution: \n");
        for (int k = 0; k < update.extent_int(0); k++) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f \n", update(k));
        }
      }
      if (linSolverStat == 1) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "NewtonFunctor: Linear solve gesv returned failure! \n");
        return NewtonSolverStatus::LinearSolveFailure;
      }

      // update solution // x = x + alpha*update
      KokkosBlas::serial_axpy(alpha, update, x);
      if (handle.debug_mode) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "NewtonFunctor: Print updated solution: \n");
        for (int k = 0; k < x.extent_int(0); k++) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f \n", x(k));
        }
      }
    }
    return NewtonSolverStatus::MaxIters;
  }  // End solve functor.
};

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // __KOKKOSBATCHED_ODE_NEWTON_HPP__
