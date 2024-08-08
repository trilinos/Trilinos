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

#ifndef KOKKOSBLAS_RUNGEKUTTATABLES_IMPL_HPP
#define KOKKOSBLAS_RUNGEKUTTATABLES_IMPL_HPP

#include <Kokkos_Array.hpp>

namespace KokkosODE {
namespace Impl {
//=====================================================================
// Generalized RK Explicit ODE solver with embedded error estimation
//=====================================================================

// Methods supported:
// Forward Euler (RKFE)
// Euler-Heun Method (RKEH)
// Fehlberg 1-2 (RKF12)
// Bogacki-Shampine (RKBS)
// Runge-Kutta 4th order (RK4)
// Fehlberg Method (RKF45)
// Cash-Karp Method (RKCK)
// Dormand-Prince Method (RKDP)

// Format follows form of Butcher Tableau
// c1| a00
// c2| a10 a11
// c3| a20 a21 a22
// c4| a30 a31 a32
// . | .   .   .
// . | .   .       .
// . | .   .          .
// cs| as0 as1 . . . . . .  ass
//--------------------------------
//   | b0  b1  b2  b3 . . . bs
//   | e0  e1  e2  e3 . . . es
//
// And is always in lower triangular form for explicit methods
// For explicit methods the methods on the diagonal will always be zero.
//
// Here, nstages = s = number of stages.
// 'order' refers to the accuracy of the method.
// The array of aij coefficients is ordered by rows as: a =
// {a00,a10,a11,a20,a21,a22....}
// e contains coefficient for error estimation

template <int order, int nstages, int variant = 0>
struct ButcherTableau {};

template <>
struct ButcherTableau<0, 0>  // Forward Euler
{
  static constexpr int order   = 1;
  static constexpr int nstages = 1;

  Kokkos::Array<double, 1> a{{1}};
  Kokkos::Array<double, nstages> b{{1}};
  Kokkos::Array<double, nstages> c{{0}};
  Kokkos::Array<double, nstages> e{{0}};
};

// Coefficients obtained from: (see page 39)
// Iserles, A.
// A First Course in the Numerical Analysis of Differential Equations."
// Cambridge: Cambridge University Press. (2008).
// https://doi:10.1017/CBO9780511995569
template <>
struct ButcherTableau<1, 1>  // Euler-Heun Method
{
  static constexpr int order   = 2;
  static constexpr int nstages = 2;  // total dimensions, nstagesxnstages system
  Kokkos::Array<double, (nstages * nstages + nstages) / 2> a{
      {0.0, 1.0, 0.0}};  //(nstages*nstages+nstages)/2 size of lower triangular matrix
  Kokkos::Array<double, nstages> b{{0.5, 0.5}};
  Kokkos::Array<double, nstages> c{{0.0, 1.0}};
  Kokkos::Array<double, nstages> e{{-0.5, 0.5}};
};

// Coefficients obtained from:
// Fehlberg, E.
// "Klassische Runge-Kutta-Formeln vierter und niedrigerer Ordnung mit
// Schrittweiten-Kontrolle und ihre Anwendung auf Wärmeleitungsprobleme."
// Computing 6, 61–71 (1970). https://doi.org/10.1007/BF02241732
template <>
struct ButcherTableau<1, 2>  // Known as Fehlberg 1-2 method
{
  static constexpr int order   = 2;
  static constexpr int nstages = 3;
  Kokkos::Array<double, (nstages * nstages + nstages) / 2> a{{0.0, 0.5, 0.0, 1.0 / 256.0, 255.0 / 256.0, 0.0}};
  Kokkos::Array<double, nstages> b{{1.0 / 512.0, 255.0 / 256.0, 1. / 512}};
  Kokkos::Array<double, nstages> c{{0.0, 1.0 / 2.0, 1.0}};
  Kokkos::Array<double, nstages> e{{1.0 / 256.0 - 1.0 / 512.0, 0.0, -1.0 / 512.0}};
};

// Coefficients obtained from:
// P. Bogacki, L.F. Shampine,
// "A 3(2) pair of Runge - Kutta formulas,"
// Applied Mathematics Letters, Volume 2, Issue 4, 1989,
// https://doi.org/10.1016/0893-9659(89)90079-7.
template <>
struct ButcherTableau<2, 3>  // Bogacki-Shampine method
{
  static constexpr int order   = 3;
  static constexpr int nstages = 4;
  Kokkos::Array<double, (nstages * nstages + nstages) / 2> a{
      {0.0, 0.5, 0.0, 0.0, 3.0 / 4.0, 0.0, 2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0}};
  Kokkos::Array<double, nstages> b{{2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0}};
  Kokkos::Array<double, nstages> c{{0.0, 0.5, 0.75, 1.0}};
  Kokkos::Array<double, nstages> e{{2.0 / 9.0 - 7.0 / 24.0, 1.0 / 3.0 - 0.25, 4.0 / 9.0 - 1.0 / 3.0, -1.0 / 8.0}};
};

// Coefficients obtained from:
// Hull, David G.
// "Fourth-order Runge-Kutta integration with stepsize control."
// AIAA Journal 15.10 (1977): 1505-1507.
template <>
struct ButcherTableau<3, 3>  // RK4
{
  static constexpr int order   = 4;
  static constexpr int nstages = 4;
  Kokkos::Array<double, (nstages * nstages + nstages) / 2> a{{0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0, 0.0}};
  Kokkos::Array<double, nstages> b{{1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0}};
  Kokkos::Array<double, nstages> c{{0.0, 0.5, 0.5, 1.0}};
  Kokkos::Array<double, nstages> e{{1.0 / 6.0, 0.0, -1.0 / 3.0, 1.0 / 6.0}};
};

// Coefficients obtained from:
// Fehlberg, E.
// "Klassische Runge-Kutta-Formeln vierter und niedrigerer Ordnung mit
// Schrittweiten-Kontrolle und ihre Anwendung auf Wärmeleitungsprobleme."
// Computing 6, 61–71 (1970). https://doi.org/10.1007/BF02241732
template <>
struct ButcherTableau<4, 5>  // Fehlberg Method
{
  static constexpr int order   = 5;
  static constexpr int nstages = 6;
  Kokkos::Array<double, (nstages * nstages + nstages) / 2> a{{0.0,
                                                              0.25,
                                                              0.0,
                                                              3.0 / 32.0,
                                                              9.0 / 32.0,
                                                              0.0,
                                                              1932.0 / 2197.0,
                                                              -7200.0 / 2197.0,
                                                              7296.0 / 2197.0,
                                                              0.0,
                                                              439.0 / 216.0,
                                                              -8.0,
                                                              3680.0 / 513.0,
                                                              -845.0 / 4104.0,
                                                              0.0,
                                                              -8.0 / 27.0,
                                                              2.0,
                                                              -3544.0 / 2565.0,
                                                              1859.0 / 4104.0,
                                                              -11.0 / 40.0,
                                                              0.0}};
  Kokkos::Array<double, nstages> b{{16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0}};
  Kokkos::Array<double, nstages> c{{0.0, 0.25, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5}};
  Kokkos::Array<double, nstages> e{{16.0 / 135.0 - 25.0 / 216.0, 0.0, 6656.0 / 12825.0 - 1408.0 / 2565.0,
                                    28561.0 / 56430.0 - 2197.0 / 4104.0, -9.0 / 50.0 + 0.2, 2.0 / 55.0}};
};

// Coefficients obtained from:
// J. R. Cash and Alan H. Karp.
// "A variable order Runge-Kutta method for initial value problems with rapidly
// varying right-hand sides." ACM Trans. Math. Softw. 16, 3 (Sept. 1990),
// 201–222. https://doi.org/10.1145/79505.79507
template <>
struct ButcherTableau<4, 5, 1>  // Cash-Karp
{
  static constexpr int order   = 5;
  static constexpr int nstages = 6;
  Kokkos::Array<double, (nstages * nstages + nstages) / 2> a{{0.0,
                                                              0.2,
                                                              0.0,
                                                              3.0 / 40.0,
                                                              9.0 / 40.0,
                                                              0.0,
                                                              0.3,
                                                              -0.9,
                                                              1.2,
                                                              0.0,
                                                              -11.0 / 54.0,
                                                              2.5,
                                                              -70.0 / 27.0,
                                                              35.0 / 27.0,
                                                              0.0,
                                                              1631.0 / 55296.0,
                                                              175.0 / 512.0,
                                                              575.0 / 13824.0,
                                                              44275.0 / 110592.0,
                                                              253.0 / 4096.0,
                                                              0.0}};
  Kokkos::Array<double, nstages> b{{37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0}};
  Kokkos::Array<double, nstages> c{{0.0, 0.2, 0.3, 0.6, 1.0, 7.0 / 8.0}};
  Kokkos::Array<double, nstages> e{{37.0 / 378.0 - 2825.0 / 27648.0, 0.0, 250.0 / 621.0 - 18575.0 / 48384.0,
                                    125.0 / 594.0 - 13525.0 / 55296.0, -277.0 / 14336.0, 512.0 / 1771.0 - 0.25}};
};

// Coefficients obtained from:
// J.R. Dormand, P.J. Prince,
// "A family of embedded Runge-Kutta formulae",
// Journal of Computational and Applied Mathematics, Volume 6, Issue 1, 1980,
// https://doi.org/10.1016/0771-050X(80)90013-3.
template <>
struct ButcherTableau<4, 6>  // Referred to as DOPRI5 or RKDP
{
  static constexpr int order   = 5;
  static constexpr int nstages = 7;
  Kokkos::Array<double, (nstages * nstages + nstages) / 2> a{{0.0,
                                                              0.2,
                                                              0.0,
                                                              3.0 / 40.0,
                                                              9.0 / 40.0,
                                                              0.0,
                                                              44.0 / 45.0,
                                                              -56.0 / 15.0,
                                                              32.0 / 9.0,
                                                              0.0,
                                                              19372.0 / 6561.0,
                                                              -25360.0 / 2187.0,
                                                              64448.0 / 6561.0,
                                                              -212.0 / 729.0,
                                                              0.0,
                                                              9017.0 / 3168.0,
                                                              -355.0 / 33.0,
                                                              46732.0 / 5247.0,
                                                              49.0 / 176.0,
                                                              -5103.0 / 18656.0,
                                                              0.0,
                                                              35.0 / 384.0,
                                                              0.0,
                                                              500.0 / 1113.0,
                                                              125.0 / 192.0,
                                                              -2187.0 / 6784.0,
                                                              11.0 / 84.0,
                                                              0.0}};
  Kokkos::Array<double, nstages> b{
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}};
  Kokkos::Array<double, nstages> c{{0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0}};
  Kokkos::Array<double, nstages> e{{35.0 / 384.0 - 5179.0 / 57600.0, 0.0, 500.0 / 1113.0 - 7571.0 / 16695.0,
                                    125.0 / 192.0 - 393.0 / 640.0, -2187.0 / 6784.0 + 92097.0 / 339200.0,
                                    11.0 / 84.0 - 187.0 / 2100.0, -1.0 / 40.0}};
};

}  // namespace Impl
}  // namespace KokkosODE

#endif  // KOKKOSBLAS_RUNGEKUTTATABLES_IMPL_HPP
