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
#ifndef __KOKKOSBATCHED_TRSM_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSM_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"
#include "KokkosBatched_InnerGemmFixA_Serial_Impl.hpp"
#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

template <typename AlgoType>
struct SerialTrsmInternalLeftLower {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
    const bool use_unit_diag, const int m, const int n, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
    const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    for (int p = 0; p < m; ++p) {
      const int iend = m - p - 1, jend = n;

      const ValueType *KOKKOS_RESTRICT a21 = iend ? A + (p + 1) * as0 + p * as1 : NULL;

      ValueType *KOKKOS_RESTRICT b1t = B + p * bs0, *KOKKOS_RESTRICT B2 = iend ? B + (p + 1) * bs0 : NULL;

      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
        for (int j = 0; j < jend; ++j) b1t[j * bs1] = b1t[j * bs1] / alpha11;
      }

      for (int i = 0; i < iend; ++i)

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
        for (int j = 0; j < jend; ++j) B2[i * bs0 + j * bs1] -= a21[i * as0] * b1t[j * bs1];
    }
  }
  return 0;
}

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(
    const bool use_unit_diag, const int m, const int n, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
    const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  constexpr int mbAlgo = Algo::Trsm::Blocked::mb();

  const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    InnerTrsmLeftLowerUnitDiag<mbAlgo> trsm_u(as0, as1, bs0, bs1);
    InnerTrsmLeftLowerNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, bs1);

    InnerGemmFixA<mbAlgo, mbAlgo> gemm(as0, as1, bs0, bs1, bs0, bs1);
    auto trsm = [&](const int ib, const int jb, const ValueType *KOKKOS_RESTRICT AA,
                    /**/ ValueType *KOKKOS_RESTRICT BB) {
      const int mb = mbAlgo;
      for (int p = 0; p < ib; p += mb) {
        const int pb = (p + mb) > ib ? (ib - p) : mb;

        // trsm update
        const ValueType *KOKKOS_RESTRICT Ap = AA + p * as0 + p * as1;
        /**/ ValueType *KOKKOS_RESTRICT Bp  = BB + p * bs0;

        if (use_unit_diag)
          trsm_u.serial_invoke(Ap, pb, jb, Bp);
        else
          trsm_n.serial_invoke(Ap, pb, jb, Bp);

        // gemm update
        for (int i = p + mb; i < ib; i += mb) {
          const int mm = (i + mb) > ib ? (ib - i) : mb;
          gemm.serial_invoke(minus_one, AA + i * as0 + p * as1, BB + p * bs0, mm, jb, pb, BB + i * bs0);
        }
      }
    };

    const bool is_small = true;  //(m*n <= 64*64);
    if (is_small) {
      trsm(m, n, A, B);
    } else {
      // // some cache blocking may need (not priority yet);
      // trsm(m, n, A, B);
    }
  }
  return 0;
}

template <typename AlgoType>
struct SerialTrsmInternalLeftUpper {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1);
};

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
    const bool use_unit_diag, const int m, const int n, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
    const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    ValueType *KOKKOS_RESTRICT B0 = B;
    for (int p = (m - 1); p >= 0; --p) {
      const int iend = p, jend = n;

      const ValueType *KOKKOS_RESTRICT a01 = A + p * as1;
      ValueType *KOKKOS_RESTRICT b1t       = B + p * bs0;

      if (!use_unit_diag) {
        const ValueType alpha11 = A[p * as0 + p * as1];

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
        for (int j = 0; j < n; ++j) b1t[j * bs1] = b1t[j * bs1] / alpha11;
      }

      if (p > 0) {  // Note: A workaround to produce correct results for
                    // complex<double> with Intel-18.2.199
        for (int i = 0; i < iend; ++i)

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
          for (int j = 0; j < jend; ++j) B0[i * bs0 + j * bs1] -= a01[i * as0] * b1t[j * bs1];
      }
    }
  }
  return 0;
}

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(
    const bool use_unit_diag, const int m, const int n, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
    const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

  constexpr int mbAlgo = Algo::Trsm::Blocked::mb();

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    InnerTrsmLeftUpperUnitDiag<mbAlgo> trsm_u(as0, as1, bs0, bs1);
    InnerTrsmLeftUpperNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, bs1);

    InnerGemmFixA<mbAlgo, mbAlgo> gemm(as0, as1, bs0, bs1, bs0, bs1);

    auto trsm = [&](const int ib, const int jb, const ValueType *KOKKOS_RESTRICT AA,
                    /**/ ValueType *KOKKOS_RESTRICT BB) {
      const int mb = mbAlgo;
      for (int pp = 0; pp < ib; pp += mb) {
        const int ptmp = ib - pp - mb, p = ptmp < 0 ? 0 : ptmp, pb = mb + (ptmp < 0) * ptmp;

        // trsm update
        const ValueType *KOKKOS_RESTRICT Ap = AA + p * as0 + p * as1;
        /**/ ValueType *KOKKOS_RESTRICT Bp  = BB + p * bs0;

        if (use_unit_diag)
          trsm_u.serial_invoke(Ap, pb, jb, Bp);
        else
          trsm_n.serial_invoke(Ap, pb, jb, Bp);

        // gemm update
        for (int i = 0; i < p; i += mb) {
          gemm.serial_invoke(minus_one, AA + i * as0 + p * as1, Bp, (i + mb) > p ? (p - i) : mb, jb, pb, BB + i * bs0);
        }
      }
    };

    const bool is_small = (m * n <= 64 * 64);
    if (is_small) {
      trsm(m, n, A, B);
    } else {
      // // some cache blocking may need (not priority yet);
      // trsm(m, n, A, B);
    }
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
