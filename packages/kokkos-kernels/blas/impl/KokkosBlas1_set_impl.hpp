// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS_SET_IMPL_HPP
#define KOKKOSBLAS_SET_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

namespace KokkosBlas {
namespace Impl {

///
/// Serial Internal Impl
/// ====================
struct SerialSetInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) A[i * as0] = alpha;

    return 0;
  }

  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (as0 > as1)
      for (int i = 0; i < m; ++i) invoke(n, alpha, A + i * as0, as1);
    else
      for (int j = 0; j < n; ++j) invoke(m, alpha, A + j * as1, as0);

    return 0;
  }
};

///
/// Team Internal Impl
/// ==================
struct TeamSetInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) { A[i * as0] = alpha; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (m > n) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m),
                           [&](const int &i) { SerialSetInternal::invoke(n, alpha, A + i * as0, as1); });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n),
                           [&](const int &j) { SerialSetInternal::invoke(m, alpha, A + j * as1, as0); });
    }
    // member.team_barrier();
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorSetInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { A[i * as0] = alpha; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (m > n) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n), [&](const int &j) { A[i * as0 + j * as1] = alpha; });
      });
    } else {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &j) { A[i * as0 + j * as1] = alpha; });
      });
    }
    // member.team_barrier();
    return 0;
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif
