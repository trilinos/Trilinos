// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_SCAL_IMPL_HPP_
#define KOKKOSBLAS1_TEAM_SCAL_IMPL_HPP_

#include <Kokkos_Core.hpp>
#include "KokkosBlas1_serial_scal_impl.hpp"

namespace KokkosBlas {
namespace Impl {

///
/// Team Internal Impl
/// ====================
struct TeamScaleInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) { A[i * as0] *= alpha; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (m > n) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m),
                           [&](const int &i) { SerialScaleInternal::invoke(n, alpha, A + i * as0, as1); });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n),
                           [&](const int &j) { SerialScaleInternal::invoke(m, alpha, A + j * as1, as0); });
    }
    // member.team_barrier();
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorScaleInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { A[i * as0] *= alpha; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (as0 > as1) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n),
                             [&](const int &j) { A[i * as0 + j * as1] *= alpha; });
      });
    } else {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &j) { A[i * as0 + j * as1] *= alpha; });
      });
    }
    // member.team_barrier();
    return 0;
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif
