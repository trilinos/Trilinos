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
#ifndef __KOKKOSBATCHED_INNER_GEMM_FIX_C_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_INNER_GEMM_FIX_C_TEAM_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_InnerGemmFixC_Decl.hpp"

namespace KokkosBatched {

template <int mb, int nb>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixC<mb, nb>::team_invoke(const MemberType &member, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int k,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, mb * nb), [&](const int &ij) {
    const int i = ij / nb, j = ij % nb;

    const ValueType *KOKKOS_RESTRICT pA = A + i * _as0, *KOKKOS_RESTRICT pB = B + j * _bs1;

    ValueType c = 0;
    for (int p = 0; p < k; ++p) c += pA[p * _as1] * pB[p * _bs0];
    C[i * _cs0 + j * _cs1] += alpha * c;
  });
  return 0;
}

template <int mb, int nb>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixC<mb, nb>::team_invoke(const MemberType &member, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int m,
                                                              const int n, const int k,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, m * n), [&](const int &ij) {
    const int i = ij / n, j = ij % n;

    const ValueType *KOKKOS_RESTRICT pA = A + i * _as0, *KOKKOS_RESTRICT pB = B + j * _bs1;

    ValueType c = 0;
    for (int p = 0; p < k; ++p) c += pA[p * _as1] * pB[p * _bs0];
    C[i * _cs0 + j * _cs1] += alpha * c;
  });
  return 0;
}
}  // namespace KokkosBatched

#endif
