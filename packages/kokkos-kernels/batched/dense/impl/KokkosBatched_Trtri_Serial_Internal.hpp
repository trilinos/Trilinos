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

#ifndef __KOKKOSBATCHED_TRTRI_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRTRI_SERIAL_INTERNAL_HPP__

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trmm_Serial_Internal.hpp"

namespace KokkosBatched {

template <typename AlgoType>
struct SerialTrtriInternalLower {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int am, const int an,
                                           ValueType *KOKKOS_RESTRICT A, const int as0, const int as1);
};

template <typename AlgoType>
struct SerialTrtriInternalUpper {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int am, const int an,
                                           ValueType *KOKKOS_RESTRICT A, const int as0, const int as1);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(const bool use_unit_diag,
                                                                                    const int am, const int /*an*/,
                                                                                    ValueType *KOKKOS_RESTRICT A,
                                                                                    const int as0, const int as1) {
  ValueType one(1.0), zero(0.0), A_ii;
  if (!use_unit_diag) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    // Check for singularity
    for (int i = 0; i < am; ++i)
      if (A[i * as0 + i * as1] == zero) return i + 1;
  }

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int i = am - 1; i >= 0; --i) {
    A[i * as0 + i * as1] = one / A[i * as0 + i * as1];

    if (i < am - 1) {
      if (use_unit_diag)
        A_ii = -one;
      else
        A_ii = -A[i * as0 + i * as1];

      ValueType *KOKKOS_RESTRICT A_subblock = &A[(i + 1) * as0 + (i + 1) * as1];
      int A_subblock_m = am - i - 1, A_subblock_n = am - i - 1;
      ValueType *KOKKOS_RESTRICT A_col_vec = &A[(i + 1) * as0 + i * as1];
      int A_col_vec_m = am - i - 1, A_col_vec_n = 1;
      // TRMV/TRMM −− x=Ax
      // A((j+1):n,j) = A((j+1):n,(j+1):n) ∗ A((j+1):n,j) ;
      SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(use_unit_diag, false, A_subblock_m, A_subblock_n,
                                                                 A_col_vec_m, A_col_vec_n, one, A_subblock, as0, as1,
                                                                 A_col_vec, as0, as1);

      // SCAL -- x=ax
      // A((j+1):n,j) = A_ii * A((j+1):n,j)
      KokkosBlas::Impl::SerialScaleInternal::invoke(A_col_vec_m, A_col_vec_n, A_ii, A_col_vec, as0, as1);
    }
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrtriInternalUpper<Algo::Trtri::Unblocked>::invoke(const bool use_unit_diag,
                                                                                    const int am, const int /*an*/,
                                                                                    ValueType *KOKKOS_RESTRICT A,
                                                                                    const int as0, const int as1) {
  ValueType one(1.0), zero(0.0), A_ii;

  if (!use_unit_diag) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    // Check for singularity
    for (int i = 0; i < am; ++i)
      if (A[i * as0 + i * as1] == zero) return i + 1;
  }

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int i = 0; i < am; ++i) {
    A[i * as0 + i * as1] = one / A[i * as0 + i * as1];

    if (i > 0) {
      if (use_unit_diag)
        A_ii = -one;
      else
        A_ii = -A[i * as0 + i * as1];

      ValueType *KOKKOS_RESTRICT A_subblock = &A[0 * as0 + 0 * as1];
      int A_subblock_m = i, A_subblock_n = i;
      ValueType *KOKKOS_RESTRICT A_col_vec = &A[0 * as0 + i * as1];
      int A_col_vec_m = i, A_col_vec_n = 1;
      // TRMV/TRMM −− x=Ax
      // A(1:(j-1),j) = A(1:(j-1),1:(j-1)) ∗ A(1:(j-1),j) ;
      // SerialTrmm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NoUnit,Algo::Trmm::Unblocked>
      SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(use_unit_diag, false, A_subblock_m, A_subblock_n,
                                                                 A_col_vec_m, A_col_vec_n, one, A_subblock, as0, as1,
                                                                 A_col_vec, as0, as1);

      // SCAL -- x=ax
      // A((j+1):n,j) = A_ii * A((j+1):n,j)
      KokkosBlas::Impl::SerialScaleInternal::invoke(A_col_vec_m, A_col_vec_n, A_ii, A_col_vec, as0, as1);
    }
  }
  return 0;
}
}  // namespace KokkosBatched
#endif  // __KOKKOSBATCHED_TRTRI_SERIAL_INTERNAL_HPP__
