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
#ifndef __KOKKOSBATCHED_SET_TRIANGULAR_INTERNAL_HPP__
#define __KOKKOSBATCHED_SET_TRIANGULAR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
struct SerialSetLowerTriangularInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const int dist, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    for (int j = 0; j < n; ++j) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = (j + dist); i < m; ++i) {
        A[i * as0 + j * as1] = alpha;
      }
    }

    return 0;
  }
};

struct TeamVectorSetLowerTriangularInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const int dist,
                                           const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &j) {
      const int jdist = j + dist;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [=](const int &i) {
        if (i >= jdist) A[i * as0 + j * as1] = alpha;
      });
    });
    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
