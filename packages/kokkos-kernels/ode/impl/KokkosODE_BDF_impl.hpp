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

#ifndef KOKKOSBLAS_BDF_IMPL_HPP
#define KOKKOSBLAS_BDF_IMPL_HPP

#include "Kokkos_Core.hpp"

#include "KokkosODE_Newton.hpp"
#include "KokkosBlas2_serial_gemv.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"

namespace KokkosODE {
namespace Impl {

template <int order>
struct BDF_table {};

template <>
struct BDF_table<1> {
  static constexpr int order = 1;
  Kokkos::Array<double, 2> coefficients{{-1.0, 1.0}};
};

template <>
struct BDF_table<2> {
  static constexpr int order = 2;
  Kokkos::Array<double, 3> coefficients{{-4.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0}};
};

template <>
struct BDF_table<3> {
  static constexpr int order = 3;
  Kokkos::Array<double, 4> coefficients{{-18.0 / 11.0, 9.0 / 11.0, -2.0 / 11.0, 6.0 / 11.0}};
};

template <>
struct BDF_table<4> {
  static constexpr int order = 4;
  Kokkos::Array<double, 5> coefficients{{-48.0 / 25.0, 36.0 / 25.0, -16.0 / 25.0, 3.0 / 25.0, 12.0 / 25.0}};
};

template <>
struct BDF_table<5> {
  static constexpr int order = 5;
  Kokkos::Array<double, 6> coefficients{
      {-300.0 / 137.0, 300.0 / 137.0, -200.0 / 137.0, 75.0 / 137.0, -12.0 / 137.0, 60.0 / 137.0}};
};

template <>
struct BDF_table<6> {
  static constexpr int order = 6;
  Kokkos::Array<double, 7> coefficients{
      {-360.0 / 147.0, 450.0 / 147.0, -400.0 / 147.0, 225.0 / 147.0, -72.0 / 147.0, 10.0 / 147.0, 60.0 / 147.0}};
};

template <class system_type, class table_type, class mv_type>
struct BDF_system_wrapper {
  const system_type mySys;
  const int neqs;
  const table_type table;
  const int order = table.order;

  double t, dt;
  mv_type yn;

  KOKKOS_FUNCTION
  BDF_system_wrapper(const system_type& mySys_, const table_type& table_, const double t_, const double dt_,
                     const mv_type& yn_)
      : mySys(mySys_), neqs(mySys_.neqs), table(table_), t(t_), dt(dt_), yn(yn_) {}

  template <class vec_type>
  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const {
    // f = f(t+dt, y)
    mySys.evaluate_function(t, dt, y, f);

    for (int eqIdx = 0; eqIdx < neqs; ++eqIdx) {
      f(eqIdx) = y(eqIdx) - table.coefficients[order] * dt * f(eqIdx);
      for (int orderIdx = 0; orderIdx < order; ++orderIdx) {
        f(eqIdx) += table.coefficients[order - 1 - orderIdx] * yn(eqIdx, orderIdx);
      }
    }
  }

  template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const {
    mySys.evaluate_jacobian(t, dt, y, jac);

    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = -table.coefficients[order] * dt * jac(rowIdx, colIdx);
      }
      jac(rowIdx, rowIdx) += 1.0;
    }
  }
};

template <class system_type, class subview_type, class d_vec_type>
struct BDF_system_wrapper2 {
  const system_type mySys;
  const int neqs;
  const subview_type psi;
  const d_vec_type d;

  bool compute_jac = true;
  double t, dt, c = 0;

  KOKKOS_FUNCTION
  BDF_system_wrapper2(const system_type& mySys_, const subview_type& psi_, const d_vec_type& d_, const double t_,
                      const double dt_)
      : mySys(mySys_), neqs(mySys_.neqs), psi(psi_), d(d_), t(t_), dt(dt_) {}

  template <class YVectorType, class FVectorType>
  KOKKOS_FUNCTION void residual(const YVectorType& y, const FVectorType& f) const {
    // f = f(t+dt, y)
    mySys.evaluate_function(t, dt, y, f);

    // std::cout << "f = psi + d - c * f = " << psi(0) << " + " << d(0) << " - "
    // << c << " * " << f(0) << std::endl;

    // rhs = higher order terms + y_{n+1}^i - y_n - dt*f
    for (int eqIdx = 0; eqIdx < neqs; ++eqIdx) {
      f(eqIdx) = psi(eqIdx) + d(eqIdx) - c * f(eqIdx);
    }
  }

  template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const {
    if (compute_jac) {
      mySys.evaluate_jacobian(t, dt, y, jac);

      // J = I - dt*(dy/dy)
      for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
        for (int colIdx = 0; colIdx < neqs; ++colIdx) {
          jac(rowIdx, colIdx) = -dt * jac(rowIdx, colIdx);
        }
        jac(rowIdx, rowIdx) += 1.0;
      }
    }
  }
};

template <class ode_type, class table_type, class vec_type, class mv_type, class mat_type, class scalar_type>
KOKKOS_FUNCTION void BDFStep(ode_type& ode, const table_type& table, scalar_type t, scalar_type dt,
                             const vec_type& y_old, const vec_type& y_new, const vec_type& rhs, const vec_type& update,
                             const vec_type& scale, const mv_type& y_vecs, const mat_type& temp, const mat_type& jac) {
  using newton_params = KokkosODE::Experimental::Newton_params;

  BDF_system_wrapper sys(ode, table, t, dt, y_vecs);
  const newton_params param(50, 1e-14, 1e-12);

  // first set y_new = y_old
  for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
    y_new(eqIdx) = y_old(eqIdx);
  }

  // solver the nonlinear problem
  { KokkosODE::Experimental::Newton::Solve(sys, param, jac, temp, y_new, rhs, update, scale); }

}  // BDFStep

template <class mat_type, class scalar_type>
KOKKOS_FUNCTION void compute_coeffs(const int order, const scalar_type factor, const mat_type& coeffs) {
  coeffs(0, 0) = 1.0;
  for (int colIdx = 0; colIdx < order; ++colIdx) {
    coeffs(0, colIdx + 1) = 1.0;
    for (int rowIdx = 0; rowIdx < order; ++rowIdx) {
      coeffs(rowIdx + 1, colIdx + 1) =
          ((rowIdx - factor * (colIdx + 1.0)) / (rowIdx + 1.0)) * coeffs(rowIdx, colIdx + 1);
    }
  }
}

template <class mat_type, class scalar_type>
KOKKOS_FUNCTION void update_D(const int order, const scalar_type factor, const mat_type& coeffs, const mat_type& tempD,
                              const mat_type& D) {
  auto subD     = Kokkos::subview(D, Kokkos::ALL(), Kokkos::pair<int, int>(0, order + 1));
  auto subTempD = Kokkos::subview(tempD, Kokkos::ALL(), Kokkos::pair<int, int>(0, order + 1));

  compute_coeffs(order, factor, coeffs);
  auto R = Kokkos::subview(coeffs, Kokkos::pair<int, int>(0, order + 1), Kokkos::pair<int, int>(0, order + 1));
  KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose,
                            KokkosBatched::Algo::Gemm::Blocked>::invoke(1.0, subD, R, 0.0, subTempD);

  compute_coeffs(order, 1.0, coeffs);
  auto U = Kokkos::subview(coeffs, Kokkos::pair<int, int>(0, order + 1), Kokkos::pair<int, int>(0, order + 1));
  KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose,
                            KokkosBatched::Algo::Gemm::Blocked>::invoke(1.0, subTempD, U, 0.0, subD);
}

template <class ode_type, class mat_type, class vec_type, class res_type, class scalar_type>
KOKKOS_FUNCTION void initial_step_size(const ode_type ode, const int order, const scalar_type t0,
                                       const scalar_type atol, const scalar_type rtol, const vec_type& y0,
                                       const res_type& f0, const mat_type& temp, scalar_type& dt_ini) {
  using KAT = Kokkos::ArithTraits<scalar_type>;

  // Extract subviews to store intermediate data
  auto scale = Kokkos::subview(temp, Kokkos::ALL(), 1);
  auto y1    = Kokkos::subview(temp, Kokkos::ALL(), 2);
  auto f1    = Kokkos::subview(temp, Kokkos::ALL(), 3);

  // Compute norms for y0 and f0
  double n0 = KAT::zero(), n1 = KAT::zero(), dt0;
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    scale(eqIdx) = atol + rtol * Kokkos::abs(y0(eqIdx));
    n0 += Kokkos::pow(y0(eqIdx) / scale(eqIdx), 2);
    n1 += Kokkos::pow(f0(eqIdx) / scale(eqIdx), 2);
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
    n2 += Kokkos::pow((f1(eqIdx) - f0(eqIdx)) / scale(eqIdx), 2);
  }
  n2 = Kokkos::sqrt(n2) / (dt0 * Kokkos::sqrt(ode.neqs));

  // Finally select initial time step dt_ini
  if ((n1 <= 1e-15) && (n2 <= 1e-15)) {
    dt_ini = Kokkos::max(1e-6, dt0 * 1e-3);
  } else {
    dt_ini = Kokkos::pow(0.01 / Kokkos::max(n1, n2), KAT::one() / (order + 1));
  }

  dt_ini = Kokkos::min(100 * dt0, dt_ini);

  // Zero out temp variables just to be safe...
  for (int eqIdx = 0; eqIdx < ode.neqs; ++eqIdx) {
    scale(eqIdx) = 0;
    y1(eqIdx)    = 0;
    f1(eqIdx)    = 0;
  }
}  // initial_step_size

template <class ode_type, class vec_type, class res_type, class mat_type, class scalar_type>
KOKKOS_FUNCTION void BDFStep(ode_type& ode, scalar_type& t, scalar_type& dt, scalar_type t_end, int& order,
                             int& num_equal_steps, const int max_newton_iters, const scalar_type atol,
                             const scalar_type rtol, const scalar_type min_factor, const vec_type& y_old,
                             const vec_type& y_new, const res_type& rhs, const res_type& update, const mat_type& temp,
                             const mat_type& temp2) {
  using newton_params = KokkosODE::Experimental::Newton_params;

  constexpr int max_order = 5;

  // For NDF coefficients see Sahmpine and Reichelt, The Matlab ODE suite, SIAM
  // SISCm 18, 1, p1-22, January 1997 Kokkos::Array<double, 6> kappa{{0.,
  // -0.1850, -1/9      , -0.0823000, -0.0415000, 0.}}; // NDF coefficients
  // kappa gamma(i) = sum_{k=1}^i(1.0 / k); gamma(0) = 0; // NDF coefficients
  // gamma_k alpha(i) = (1 - kappa(i)) * gamma(i) error_const(i) = kappa(i) *
  // gamma(i) + 1 / (i + 1)
  const Kokkos::Array<const double, 6> alpha{{0., 1.185, 1.66666667, 1.98421667, 2.16979167, 2.28333333}};
  const Kokkos::Array<const double, 6> error_const{{1., 0.315, 0.16666667, 0.09911667, 0.11354167, 0.16666667}};

  // Extract columns of temp to form temporary
  // subviews to operate on.
  // const int numRows = temp.extent_int(0); const int numCols =
  // temp.extent_int(1); std::cout << "numRows: " << numRows << ", numCols: " <<
  // numCols << std::endl; std::cout << "Extract subview from temp" <<
  // std::endl;
  int offset = 2;
  auto D = Kokkos::subview(temp, Kokkos::ALL(), Kokkos::pair<int, int>(offset, offset + 8));  // y and its derivatives
  offset += 8;
  auto tempD = Kokkos::subview(temp, Kokkos::ALL(), Kokkos::pair<int, int>(offset, offset + 8));
  offset += 8;
  auto scale = Kokkos::subview(temp, Kokkos::ALL(), offset + 1);
  ++offset;  // Scaling coefficients for error calculation
  auto y_predict = Kokkos::subview(temp, Kokkos::ALL(), offset + 1);
  ++offset;  // Initial guess for y_{n+1}
  auto psi = Kokkos::subview(temp, Kokkos::ALL(), offset + 1);
  ++offset;  // Higher order terms contribution to rhs
  auto error = Kokkos::subview(temp, Kokkos::ALL(), offset + 1);
  ++offset;  // Error estimate
  auto jac =
      Kokkos::subview(temp, Kokkos::ALL(), Kokkos::pair<int, int>(offset, offset + ode.neqs));  // Jacobian matrix
  offset += ode.neqs;
  auto tmp_gesv = Kokkos::subview(
      temp, Kokkos::ALL(), Kokkos::pair<int, int>(offset, offset + ode.neqs + 4));  // Buffer space for gesv calculation
  offset += ode.neqs + 4;

  auto coeffs = Kokkos::subview(temp2, Kokkos::ALL(), Kokkos::pair<int, int>(0, 6));
  auto gamma  = Kokkos::subview(temp2, Kokkos::ALL(), 6);
  gamma(0)    = 0.0;
  gamma(1)    = 1.0;
  gamma(2)    = 1.5;
  gamma(3)    = 1.83333333;
  gamma(4)    = 2.08333333;
  gamma(5)    = 2.28333333;

  BDF_system_wrapper2 sys(ode, psi, update, t, dt);
  const newton_params param(
      max_newton_iters, atol,
      Kokkos::max(10 * Kokkos::ArithTraits<scalar_type>::eps() / rtol, Kokkos::min(0.03, Kokkos::sqrt(rtol))));

  scalar_type max_step = Kokkos::ArithTraits<scalar_type>::max();
  scalar_type min_step = Kokkos::ArithTraits<scalar_type>::min();
  scalar_type safety   = 0.675, error_norm;
  if (dt > max_step) {
    update_D(order, max_step / dt, coeffs, tempD, D);
    dt              = max_step;
    num_equal_steps = 0;
  } else if (dt < min_step) {
    update_D(order, min_step / dt, coeffs, tempD, D);
    dt              = min_step;
    num_equal_steps = 0;
  }

  // first set y_new = y_old
  for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
    y_new(eqIdx) = y_old(eqIdx);
  }

  double t_new       = 0;
  bool step_accepted = false;
  while (!step_accepted) {
    if (dt < min_step) {
      return;
    }
    t_new = t + dt;

    if (t_new > t_end) {
      t_new = t_end;
      update_D(order, (t_new - t) / dt, coeffs, tempD, D);
      num_equal_steps = 0;
    }
    dt = t_new - t;

    for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
      y_predict(eqIdx) = 0;
      for (int orderIdx = 0; orderIdx < order + 1; ++orderIdx) {
        y_predict(eqIdx) += D(eqIdx, orderIdx);
      }
      scale(eqIdx) = atol + rtol * Kokkos::abs(y_predict(eqIdx));
    }

    // Compute psi, the sum of the higher order
    // contribution to the residual
    auto subD     = Kokkos::subview(D, Kokkos::ALL(), Kokkos::pair<int, int>(1, order + 1));
    auto subGamma = Kokkos::subview(gamma, Kokkos::pair<int, int>(1, order + 1));
    KokkosBlas::Experimental::serial_gemv('N', 1.0 / alpha[order], subD, subGamma, 0.0, psi);

    sys.compute_jac = true;
    sys.c           = dt / alpha[order];
    sys.jacobian(y_new, jac);
    sys.compute_jac = true;
    Kokkos::Experimental::local_deep_copy(y_new, y_predict);
    Kokkos::Experimental::local_deep_copy(update, 0);
    KokkosODE::Experimental::newton_solver_status newton_status =
        KokkosODE::Experimental::Newton::Solve(sys, param, jac, tmp_gesv, y_new, rhs, update, scale);

    for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
      update(eqIdx) = y_new(eqIdx) - y_predict(eqIdx);
    }

    if (newton_status == KokkosODE::Experimental::newton_solver_status::MAX_ITER) {
      dt = 0.5 * dt;
      update_D(order, 0.5, coeffs, tempD, D);
      num_equal_steps = 0;

    } else {
      // Estimate the solution error
      safety     = 0.9 * (2 * max_newton_iters + 1) / (2 * max_newton_iters + param.iters);
      error_norm = 0;
      for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
        scale(eqIdx) = atol + rtol * Kokkos::abs(y_new(eqIdx));
        error(eqIdx) = error_const[order] * update(eqIdx) / scale(eqIdx);
        error_norm += error(eqIdx) * error(eqIdx);
      }
      error_norm = Kokkos::sqrt(error_norm) / Kokkos::sqrt(sys.neqs);

      // Check error norm and adapt step size or accept step
      if (error_norm > 1) {
        scalar_type factor = Kokkos::max(min_factor, safety * Kokkos::pow(error_norm, -1.0 / (order + 1)));
        dt                 = factor * dt;
        update_D(order, factor, coeffs, tempD, D);
        num_equal_steps = 0;
      } else {
        step_accepted = true;
      }
    }
  }  // while(!step_accepted)

  // Now that our time step has been
  // accepted we update all our states
  // and see if we can adapt the order
  // or the time step before going to
  // the next step.
  ++num_equal_steps;
  t = t_new;
  for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
    D(eqIdx, order + 2) = update(eqIdx) - D(eqIdx, order + 1);
    D(eqIdx, order + 1) = update(eqIdx);
    for (int orderIdx = order; 0 <= orderIdx; --orderIdx) {
      D(eqIdx, orderIdx) += D(eqIdx, orderIdx + 1);
    }
  }

  // Not enough steps at constant dt
  // have been succesfull so we do not
  // attempt order adaptation.
  double error_low = 0, error_high = 0;
  if (num_equal_steps < order + 1) {
    return;
  }

  if (1 < order) {
    for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
      error_low += Kokkos::pow(error_const[order - 1] * D(eqIdx, order) / scale(eqIdx), 2);
    }
    error_low = Kokkos::sqrt(error_low) / Kokkos::sqrt(sys.neqs);
  } else {
    error_low = Kokkos::ArithTraits<double>::max();
  }

  if (order < max_order) {
    for (int eqIdx = 0; eqIdx < sys.neqs; ++eqIdx) {
      error_high += Kokkos::pow(error_const[order + 1] * D(eqIdx, order + 2) / scale(eqIdx), 2);
    }
    error_high = Kokkos::sqrt(error_high) / Kokkos::sqrt(sys.neqs);
  } else {
    error_high = Kokkos::ArithTraits<double>::max();
  }

  double factor_low, factor_mid, factor_high, factor;
  factor_low  = Kokkos::pow(error_low, -1.0 / order);
  factor_mid  = Kokkos::pow(error_norm, -1.0 / (order + 1));
  factor_high = Kokkos::pow(error_high, -1.0 / (order + 2));

  int delta_order = 0;
  if ((factor_mid < factor_low) && (factor_high < factor_low)) {
    delta_order = -1;
    factor      = factor_low;
  } else if ((factor_low < factor_high) && (factor_mid < factor_high)) {
    delta_order = 1;
    factor      = factor_high;
  } else {
    delta_order = 0;
    factor      = factor_mid;
  }
  order += delta_order;
  factor = Kokkos::fmin(10, safety * factor);
  dt *= factor;

  update_D(order, factor, coeffs, tempD, D);
  num_equal_steps = 0;

}  // BDFStep

}  // namespace Impl
}  // namespace KokkosODE

#endif  // KOKKOSBLAS_BDF_IMPL_HPP
