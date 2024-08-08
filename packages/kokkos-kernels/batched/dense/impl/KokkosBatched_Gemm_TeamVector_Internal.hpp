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
#ifndef __KOKKOSBATCHED_GEMM_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMM_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_team_scal_impl.hpp"

namespace KokkosBatched {

///
/// TeamVector Internal Impl
/// ====================
template <typename ArgAlgo, bool useConjA = false>
struct TeamVectorGemmInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const int k,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                           const int bs1, const ScalarType beta,
                                           /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1);
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorGemmInternal<Algo::Gemm::Unblocked, false>::invoke(
    const MemberType &member, const int m, const int n, const int k, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
    const int bs1, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1) {
  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)

  const ScalarType one(1.0), zero(0.0);

  if (beta == zero)
    KokkosBlas::Impl::TeamVectorSetInternal::invoke(member, m, n, zero, C, cs0, cs1);
  else if (beta != one)
    KokkosBlas::Impl::TeamVectorScaleInternal::invoke(member, m, n, beta, C, cs0, cs1);

  if (alpha != ScalarType(0.0)) {
    if (m <= 0 || n <= 0 || k <= 0) return 0;

    if (beta != one) member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
      const ValueType *KOKKOS_RESTRICT pA = A + i * as0;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n), [&](const int &j) {
        const ValueType *KOKKOS_RESTRICT pB = B + j * bs1;

        ValueType c = ValueType(0);
        for (int p = 0; p < k; ++p) c += pA[p * as1] * pB[p * bs0];
        C[i * cs0 + j * cs1] += alpha * c;
      });
    });
  }
  return 0;
}

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorGemmInternal<Algo::Gemm::Unblocked, true>::invoke(
    const MemberType &member, const int m, const int n, const int k, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
    const int bs1, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1) {
  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)

  const ScalarType one(1.0), zero(0.0);

  if (beta == zero)
    KokkosBlas::Impl::TeamVectorSetInternal::invoke(member, m, n, zero, C, cs0, cs1);
  else if (beta != one)
    KokkosBlas::Impl::TeamVectorScaleInternal::invoke(member, m, n, beta, C, cs0, cs1);

  if (alpha != ScalarType(0.0)) {
    if (m <= 0 || n <= 0 || k <= 0) return 0;

    if (beta != one) member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
      const ValueType *KOKKOS_RESTRICT pA = A + i * as0;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n), [&](const int &j) {
        const ValueType *KOKKOS_RESTRICT pB = B + j * bs1;

        ValueType c = ValueType(0);
        for (int p = 0; p < k; ++p) c += Kokkos::ArithTraits<ValueType>::conj(pA[p * as1]) * pB[p * bs0];
        C[i * cs0 + j * cs1] += alpha * c;
      });
    });
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
