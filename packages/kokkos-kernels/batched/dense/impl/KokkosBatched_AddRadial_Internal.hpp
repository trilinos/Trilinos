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
#ifndef __KOKKOSBATCHED_ADD_RADIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_ADD_RADIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
struct SerialAddRadialInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ScalarType tiny,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as) {
    const auto abs_tiny       = tiny > 0 ? tiny : -tiny;
    const auto minus_abs_tiny = -abs_tiny;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      // const auto a_real = RealPart(A[i*as]);
      const auto a_real = Kokkos::ArithTraits<ValueType>::real(A[i * as]);
      A[i * as] += ValueType(minus_abs_tiny) * ValueType(a_real < 0);
      A[i * as] += ValueType(abs_tiny) * ValueType(a_real >= 0);
    }

    return 0;
  }
};

///
/// Team Internal Impl
/// ==================
struct TeamAddRadialInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const ScalarType tiny,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as) {
    const auto abs_tiny       = tiny > 0 ? tiny : -tiny;
    const auto minus_abs_tiny = -abs_tiny;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
      // const auto a_real = RealPart(A[i*as]);
      const auto a_real = Kokkos::ArithTraits<ValueType>::real(A[i * as]);
      A[i * as] += ValueType(minus_abs_tiny) * ValueType(a_real < 0);
      A[i * as] += ValueType(abs_tiny) * ValueType(a_real >= 0);
    });

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
