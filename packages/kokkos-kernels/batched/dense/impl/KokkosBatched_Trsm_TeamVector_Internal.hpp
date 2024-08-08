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
#ifndef __KOKKOSBATCHED_TRSM_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSM_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_team_scal_impl.hpp"

namespace KokkosBatched {

///
/// Team Internal Impl
/// ====================

template <typename AlgoType>
struct TeamVectorTrsmInternalLeftLower {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const bool use_unit_diag, const int m, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::TeamVectorSetInternal::invoke(member, m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamVectorScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    for (int p = 0; p < m; ++p) {
      // Made this non-const in order to WORKAROUND issue #349
      int iend = m - p - 1;
      int jend = n;

      const ValueType *KOKKOS_RESTRICT a21 = iend ? A + (p + 1) * as0 + p * as1 : NULL;

      ValueType *KOKKOS_RESTRICT b1t = B + p * bs0, *KOKKOS_RESTRICT B2 = iend ? B + (p + 1) * bs0 : NULL;

      member.team_barrier();
      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, jend),
                             [&](const int &j) { b1t[j * bs1] = b1t[j * bs1] / alpha11; });
        member.team_barrier();
      }
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, iend), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, jend), [&](const int &j) {
          // assume layout right for batched computation
          B2[i * bs0 + j * bs1] -= a21[i * as0] * b1t[j * bs1];
        });
      });
    }
  }
  return 0;
}

template <typename AlgoType>
struct TeamVectorTrsmInternalLeftUpper {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const bool use_unit_diag, const int m, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
    const MemberType &member, const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);

  // note that parallel range is different ( m*n vs m-1*n);
  if (alpha == zero)
    KokkosBlas::Impl::TeamVectorSetInternal::invoke(member, m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::TeamVectorScaleInternal::invoke(member, m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    ValueType *KOKKOS_RESTRICT B0 = B;
    for (int p = (m - 1); p >= 0; --p) {
      // Made this non-const in order to WORKAROUND issue #349
      int iend = p;
      int jend = n;

      const ValueType *KOKKOS_RESTRICT a01 = A + p * as1;
      /**/ ValueType *KOKKOS_RESTRICT b1t  = B + p * bs0;

      member.team_barrier();
      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, jend),
                             [&](const int &j) { b1t[j * bs1] = b1t[j * bs1] / alpha11; });
        member.team_barrier();
      }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, iend), [&](const int &i) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, jend),
                             [&](const int &j) { B0[i * bs0 + j * bs1] -= a01[i * as0] * b1t[j * bs1]; });
      });
    }
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
