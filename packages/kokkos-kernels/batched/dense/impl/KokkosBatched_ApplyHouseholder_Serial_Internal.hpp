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
#ifndef __KOKKOSBATCHED_APPLY_HOUSEHOLDER_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_APPLY_HOUSEHOLDER_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialApplyLeftHouseholderInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ValueType* tau,
                                           /* */ ValueType* u2, const int u2s,
                                           /* */ ValueType* a1t, const int a1ts,
                                           /* */ ValueType* A2, const int as0, const int as1,
                                           /* */ ValueType* w1t) {
    typedef ValueType value_type;

    /// u2  m x 1
    /// a1t 1 x n
    /// A2  m x n

    // apply a single householder transform H from the left to a row vector a1t
    // and a matrix A2
    const value_type inv_tau = value_type(1) / (*tau);

    // compute the followings:
    // a1t -=    inv(tau)(a1t + u2'A2)
    // A2  -= u2 inv(tau)(a1t + u2'A2)

    // w1t = a1t + u2'A2 = A2^T conj(u2)
    // w1t /= tau
    for (int j = 0; j < n; ++j) {
      value_type tmp = a1t[j * a1ts];
      for (int i = 0; i < m; ++i) tmp += Kokkos::ArithTraits<value_type>::conj(u2[i * u2s]) * A2[i * as0 + j * as1];
      w1t[j] = tmp * inv_tau;  // /= (*tau);
    }

    // a1t -= w1t    (axpy)
    for (int j = 0; j < n; ++j) a1t[j * a1ts] -= w1t[j];

    // A2  -= u2 w1t (ger)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i) A2[i * as0 + j * as1] -= u2[i * u2s] * w1t[j];

    return 0;
  }
};

struct SerialApplyRightHouseholderInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ValueType* tau,
                                           /* */ ValueType* u2, const int u2s,
                                           /* */ ValueType* a1, const int a1s,
                                           /* */ ValueType* A2, const int as0, const int as1,
                                           /* */ ValueType* w1) {
    typedef ValueType value_type;
    /// u2 n x 1
    /// a1 m x 1
    /// A2 m x n

    // apply a single householder transform H from the left to a row vector a1t
    // and a matrix A2
    const value_type inv_tau = value_type(1) / (*tau);

    // compute the followings:
    // a1 -= inv(tau)(a1 + A2 u2)
    // A2 -= inv(tau)(a1 + A2 u2) u2'

    // w1 = a1 + A2 u2
    // w1 /= tau
    for (int i = 0; i < m; ++i) {
      value_type tmp = a1[i * a1s];
      for (int j = 0; j < n; ++j) tmp += A2[i * as0 + j * as1] * u2[j * u2s];
      w1[i] = tmp * inv_tau;  // \= (*tau);
    }

    // a1 -= w1 (axpy)
    for (int i = 0; i < m; ++i) a1[i * a1s] -= w1[i];

    // A2 -= w1 * u2' (ger with conjugate)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i) A2[i * as0 + j * as1] -= w1[i] * Kokkos::ArithTraits<ValueType>::conj(u2[j * u2s]);

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
