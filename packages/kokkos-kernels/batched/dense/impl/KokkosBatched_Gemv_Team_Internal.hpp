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
#ifndef __KOKKOSBATCHED_GEMV_TEAM_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMV_TEAM_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

// #include "KokkosBlas1_set_impl.hpp"
// #include "KokkosBlas1_team_scal_impl.hpp"
// #include "KokkosBlas2_serial_gemv_inner_multiple_dot.hpp"

namespace KokkosBatched {

///
/// Team Internal Impl
/// ====================
template <typename ArgAlgo>
struct TeamGemvInternal {
  template <typename MemberType, typename ScalarType, typename layout, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int N, const int m, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1, const int as2, const ValueType *KOKKOS_RESTRICT x,
                                           const int xs0, const int xs1, const ScalarType beta,
                                           /**/ ValueType *KOKKOS_RESTRICT y, const int ys0, const int ys1);
};

template <>
template <typename MemberType, typename ScalarType, typename layout, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamGemvInternal<Algo::Gemv::Unblocked>::invoke(
    const MemberType &member, const int N, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1, const int as2, const ValueType *KOKKOS_RESTRICT X,
    const int xs0, const int xs1, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT Y, const int ys0, const int ys1) {
  const ScalarType one(1.0), zero(0.0);

  // y_l = beta y_l + alpha A_l x_l for l in range(0, N)
  // y_l (m), A_l(m x n), B_l(n)

  if (beta == zero)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, N * m), [&](const int &iTemp) {
      int iRow, iMatrix;
      getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
      Y[ys0 * iMatrix + ys1 * iRow] = zero;
    });
  else if (beta != one)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, N * m), [&](const int &iTemp) {
      int iRow, iMatrix;
      getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
      Y[ys0 * iMatrix + ys1 * iRow] *= beta;
    });

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    if (beta != one) member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, N * m), [&](const int &iTemp) {
      int iRow, iMatrix;
      ValueType t(0);
      getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
      for (int i = 0; i < n; ++i) t += A[as0 * iMatrix + as1 * iRow + as2 * i] * X[xs0 * iMatrix + xs1 * i];
      Y[ys0 * iMatrix + ys1 * iRow] += alpha * t;
    });
  }
  return 0;
}
}  // namespace KokkosBatched

#endif
