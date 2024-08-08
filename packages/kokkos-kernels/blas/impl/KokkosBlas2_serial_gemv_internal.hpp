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
#ifndef __KOKKOSBLAS_GEMV_SERIAL_INTERNAL_HPP__
#define __KOKKOSBLAS_GEMV_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBlas_util.hpp"
#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"
#include "KokkosBlas2_serial_gemv_inner_multiple_dot.hpp"

namespace KokkosBlas {
namespace Impl {
///
/// Serial Internal Impl
/// ====================

template <typename ArgAlgo>
struct SerialGemvInternal {
  template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
  KOKKOS_INLINE_FUNCTION static int invoke(OpA op, const int m, const int n, const ScalarType alpha,
                                           const ValueAType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           const ValueXType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
                                           /**/ ValueYType *KOKKOS_RESTRICT y, const int ys0);

  // default OpA = OpID
  template <typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ScalarType alpha,
                                           const ValueAType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           const ValueXType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
                                           /**/ ValueYType *KOKKOS_RESTRICT y, const int ys0) {
    return invoke(OpID(), m, n, alpha, A, as0, as1, x, xs0, beta, y, ys0);
  }
};

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int SerialGemvInternal<Algo::Gemv::Unblocked>::invoke(
    OpA op, const int m, const int n, const ScalarType alpha, const ValueAType *KOKKOS_RESTRICT A, const int as0,
    const int as1, const ValueXType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
    /**/ ValueYType *KOKKOS_RESTRICT y, const int ys0) {
  const ScalarType one(1.0), zero(0.0);

  // y = beta y + alpha A x
  // y (m), A(m x n), B(n)

  if (beta == zero) {
    ValueYType val_zero{};  // can be Vector<SIMD> so avoid assigning explicit 0
    KokkosBlas::Impl::SerialSetInternal::invoke(m, val_zero, y, ys0);
  } else if (beta != one)
    KokkosBlas::Impl::SerialScaleInternal::invoke(m, beta, y, ys0);

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    for (int i = 0; i < m; ++i) {
      ValueYType t(0);
      const ValueAType *KOKKOS_RESTRICT tA = A + i * as0;

#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int j = 0; j < n; ++j) t += op(tA[j * as1]) * x[j * xs0];
      y[i * ys0] += alpha * t;
    }
  }
  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int SerialGemvInternal<Algo::Gemv::Blocked>::invoke(
    OpA /* op */, const int m, const int n, const ScalarType alpha, const ValueAType *KOKKOS_RESTRICT A, const int as0,
    const int as1, const ValueXType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
    /**/ ValueYType *KOKKOS_RESTRICT y, const int ys0) {
  const ScalarType one(1.0), zero(0.0);

  // y = beta y + alpha A x
  // y (m), A(m x n), B(n)

  constexpr int mbAlgo = Algo::Gemv::Blocked::mb();

  if (beta == zero)
    Impl::SerialSetInternal::invoke(m, zero, y, ys0);
  else if (beta != one)
    Impl::SerialScaleInternal::invoke(m, beta, y, ys0);

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    Impl::InnerMultipleDotProduct<mbAlgo> inner(as0, as1, xs0, ys0);
    const int mb = mbAlgo;
    for (int i = 0; i < m; i += mb)
      inner.serial_invoke<OpA>(alpha, A + i * as0, x, (i + mb) > m ? (m - i) : mb, n, y + i * ys0);
  }
  return 0;
}
}  // namespace Impl
}  // namespace KokkosBlas

#endif
