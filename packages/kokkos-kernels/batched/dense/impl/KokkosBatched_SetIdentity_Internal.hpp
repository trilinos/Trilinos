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
#ifndef __KOKKOSBATCHED_SET_IDENTITY_INTERNAL_HPP__
#define __KOKKOSBATCHED_SET_IDENTITY_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
struct SerialSetIdentityInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    const ValueType one(1), zero(0);
    for (int j = 0; j < n; ++j) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = 0; i < m; ++i) {
        A[i * as0 + j * as1] = i == j ? one : zero;
      }
    }

    return 0;
  }
};

///
/// Team Internal Impl
/// ==================
struct TeamSetIdentityInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    const ValueType one(1), zero(0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int j = 0; j < n; ++j) A[i * as0 + j * as1] = i == j ? one : zero;
    });

    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorSetIdentityInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    const ValueType one(1), zero(0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n),
                           [&](const int &j) { A[i * as0 + j * as1] = i == j ? one : zero; });
    });

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
