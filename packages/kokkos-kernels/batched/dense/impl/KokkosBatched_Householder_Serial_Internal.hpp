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
#ifndef __KOKKOSBATCHED_HOUSEHOLDER_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_HOUSEHOLDER_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialLeftHouseholderInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m_x2,
                                           /* */ ValueType* chi1,
                                           /* */ ValueType* x2, const int x2s,
                                           /* */ ValueType* tau) {
    typedef ValueType value_type;
    typedef typename Kokkos::ArithTraits<ValueType>::mag_type mag_type;

    const mag_type zero(0);
    const mag_type half(0.5);
    const mag_type one(1);
    const mag_type minus_one(-1);

    /// compute the 2norm of x2
    mag_type norm_x2_square(0);
    for (int i = 0; i < m_x2; ++i) {
      const auto x2_at_i = x2[i * x2s];
      norm_x2_square += x2_at_i * x2_at_i;
    }

    /// if norm_x2 is zero, return with trivial values
    if (norm_x2_square == zero) {
      *chi1 = -(*chi1);
      *tau  = half;

      return 0;
    }

    /// compute magnitude of chi1, equal to norm2 of chi1
    const mag_type norm_chi1 = Kokkos::ArithTraits<value_type>::abs(*chi1);

    /// compute 2 norm of x using norm_chi1 and norm_x2
    const mag_type norm_x = Kokkos::ArithTraits<mag_type>::sqrt(norm_x2_square + norm_chi1 * norm_chi1);

    /// compute alpha
    const mag_type alpha = (*chi1 < 0 ? one : minus_one) * norm_x;

    /// overwrite x2 with u2
    const value_type chi1_minus_alpha     = *chi1 - alpha;
    const value_type inv_chi1_minus_alpha = one / chi1_minus_alpha;
    for (int i = 0; i < m_x2; ++i) x2[i * x2s] *= inv_chi1_minus_alpha;

    // later consider to use the following
    // SerialScaleInternal::invoke(m_x2, inv_chi1_minus_alpha, x2, x2s);

    /// compute tau
    const mag_type chi1_minus_alpha_square = chi1_minus_alpha * chi1_minus_alpha;
    *tau                                   = half + half * (norm_x2_square / chi1_minus_alpha_square);

    /// overwrite chi1 with alpha
    *chi1 = alpha;

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
