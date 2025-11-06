// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_COPY_INTERNAL_HPP
#define KOKKOSBATCHED_COPY_INTERNAL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

struct SerialCopyInternal {
  template <typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const int m, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                /* */ ValueType *KOKKOS_RESTRICT B, const int bs0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) B[i * bs0] = A[i * as0];

    return 0;
  }
  template <typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const int m, const int n, const ValueType *KOKKOS_RESTRICT A,
                                                const int as0, const int as1,
                                                /* */ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
    if (as1 < as0)
      for (int i = 0; i < m; ++i) invoke(n, A + i * as0, as1, B + i * bs0, bs1);
    else
      for (int j = 0; j < n; ++j) invoke(m, A + j * as1, as0, B + j * bs1, bs0);
    return 0;
  }
};

///
/// Team Internal Impl
/// ==================
struct TeamCopyInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                /* */ ValueType *KOKKOS_RESTRICT B, const int bs0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) { B[i * bs0] = A[i * as0]; });
    // member.team_barrier();
    return 0;
  }
  template <typename MemberType, typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                /* */ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
    if (m >= n) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m),
                           [&](const int &i) { SerialCopyInternal::invoke(n, A + i * as0, as1, B + i * bs0, bs1); });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n),
                           [&](const int &j) { SerialCopyInternal::invoke(m, A + j * as1, as0, B + j * bs1, bs0); });
    }
    // member.team_barrier();
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorCopyInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                /* */ ValueType *KOKKOS_RESTRICT B, const int bs0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { B[i * bs0] = A[i * as0]; });
    // member.team_barrier();
    return 0;
  }
  template <typename MemberType, typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                /* */ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
    if (as0 > as1) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n),
                             [&](const int &j) { B[i * bs0 + j * bs1] = A[i * as0 + j * as1]; });
      });
    } else {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n),
                             [&](const int &j) { B[i * bs0 + j * bs1] = A[i * as0 + j * as1]; });
      });
    }
    // member.team_barrier();
    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
