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
#ifndef __KOKKOSBATCHED_TRSV_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSV_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"
#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"
#include "KokkosBlas2_serial_gemv_internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

///
/// Lower
///

template <typename AlgoType>
struct SerialTrsvInternalLower {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int m, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT b, const int bs0);
};

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(const bool use_unit_diag, const int m,
                                                                                  const ScalarType alpha,
                                                                                  const ValueType *KOKKOS_RESTRICT A,
                                                                                  const int as0, const int as1,
                                                                                  /**/ ValueType *KOKKOS_RESTRICT b,
                                                                                  const int bs0) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, zero, b, bs0);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, alpha, b, bs0);
    if (m <= 0) return 0;

    for (int p = 0; p < m; ++p) {
      const int iend = m - p - 1;

      const ValueType *KOKKOS_RESTRICT a21 = iend ? A + (p + 1) * as0 + p * as1 : NULL;

      ValueType *KOKKOS_RESTRICT beta1 = b + p * bs0, *KOKKOS_RESTRICT b2 = iend ? beta1 + bs0 : NULL;

      // with KOKKOS_RESTRICT a compiler assumes that the pointer is not
      // accessed by others op(/=) uses this pointer and changes the associated
      // values, which brings a compiler problem
      if (!use_unit_diag) *beta1 = *beta1 / A[p * as0 + p * as1];

      for (int i = 0; i < iend; ++i) b2[i * bs0] -= a21[i * as0] * (*beta1);
    }
  }
  return 0;
}

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsvInternalLower<Algo::Trsv::Blocked>::invoke(const bool use_unit_diag, const int m,
                                                                                const ScalarType alpha,
                                                                                const ValueType *KOKKOS_RESTRICT A,
                                                                                const int as0, const int as1,
                                                                                /**/ ValueType *KOKKOS_RESTRICT b,
                                                                                const int bs0) {
  const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

  constexpr int mbAlgo = Algo::Trsv::Blocked::mb();

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, zero, b, bs0);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, alpha, b, bs0);
    if (m <= 0) return 0;

    /// case GPU: team size is large and blocksize (mb,nb) is small
    InnerTrsmLeftLowerUnitDiag<mbAlgo> trsm_u(as0, as1, bs0, 0);
    InnerTrsmLeftLowerNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, 0);

    const int mb = mbAlgo;
    for (int p = 0; p < m; p += mb) {
      const int pb = ((p + mb) > m ? (m - p) : mb);

      // trsm update
      const ValueType *KOKKOS_RESTRICT Ap = A + p * as0 + p * as1;
      /**/ ValueType *KOKKOS_RESTRICT bp  = b + p * bs0;

      if (use_unit_diag)
        trsm_u.serial_invoke(Ap, pb, 1, bp);
      else
        trsm_n.serial_invoke(Ap, pb, 1, bp);

      // gemv update
      KokkosBlas::Impl::SerialGemvInternal<Algo::Gemv::Blocked>::invoke(m - p - pb, pb, minus_one, Ap + pb * as0, as0,
                                                                        as1, bp, bs0, one, bp + pb * bs0, bs0);
    }
  }
  return 0;
}

///
/// Upper
///

template <typename AlgoType>
struct SerialTrsvInternalUpper {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int m, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT b, const int bs0);
};

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(const bool use_unit_diag, const int m,
                                                                                  const ScalarType alpha,
                                                                                  const ValueType *KOKKOS_RESTRICT A,
                                                                                  const int as0, const int as1,
                                                                                  /**/ ValueType *KOKKOS_RESTRICT b,
                                                                                  const int bs0) {
  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, zero, b, bs0);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, alpha, b, bs0);
    if (m <= 0) return 0;

    ValueType *KOKKOS_RESTRICT b0 = b;
    for (int p = (m - 1); p >= 0; --p) {
      const int iend = p;

      const ValueType *KOKKOS_RESTRICT a01  = A + p * as1;
      /**/ ValueType *KOKKOS_RESTRICT beta1 = b + p * bs0;

      // with KOKKOS_RESTRICT a compiler assumes that the pointer is not
      // accessed by others op(/=) uses this pointer and changes the associated
      // values, which brings a compiler problem
      if (!use_unit_diag) *beta1 = *beta1 / A[p * as0 + p * as1];

      for (int i = 0; i < iend; ++i) b0[i * bs0] -= a01[i * as0] * (*beta1);
    }
  }
  return 0;
}

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(const bool use_unit_diag, const int m,
                                                                                const ScalarType alpha,
                                                                                const ValueType *KOKKOS_RESTRICT A,
                                                                                const int as0, const int as1,
                                                                                /**/ ValueType *KOKKOS_RESTRICT b,
                                                                                const int bs0) {
  const ScalarType one(1.0), zero(0.0), minus_one(-1.0);

  constexpr int mbAlgo = Algo::Trsm::Blocked::mb();

  // note that parallel range is different ( m*n vs m-1*n);
  if (alpha == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, zero, b, bs0);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, alpha, b, bs0);
    if (m <= 0) return 0;

    InnerTrsmLeftUpperUnitDiag<mbAlgo> trsm_u(as0, as1, bs0, 0);
    InnerTrsmLeftUpperNonUnitDiag<mbAlgo> trsm_n(as0, as1, bs0, 0);

    const int mb = mbAlgo;
    for (int pp = 0; pp < m; pp += mb) {
      const int ptmp = (m - pp - mb), p = (ptmp < 0 ? 0 : ptmp), pb = (mb + (ptmp < 0) * ptmp);

      // trsm update
      const ValueType *KOKKOS_RESTRICT Ap = A + p * as0 + p * as1;
      /**/ ValueType *KOKKOS_RESTRICT bp  = b + p * bs0;

      if (use_unit_diag)
        trsm_u.serial_invoke(Ap, pb, 1, bp);
      else
        trsm_n.serial_invoke(Ap, pb, 1, bp);

      // gemv update
      KokkosBlas::Impl::SerialGemvInternal<Algo::Gemv::Blocked>::invoke(p, pb, minus_one, Ap - p * as0, as0, as1, bp,
                                                                        bs0, one, b, bs0);
    }
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
