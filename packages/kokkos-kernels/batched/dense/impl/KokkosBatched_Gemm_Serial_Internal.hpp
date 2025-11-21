// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMM_SERIAL_INTERNAL_HPP
#define KOKKOSBATCHED_GEMM_SERIAL_INTERNAL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"

#include "KokkosBatched_InnerGemmFixC_Serial_Impl.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

template <typename ArgAlgo>
struct SerialGemmInternal {
  template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(OpA opA, OpB opB, const int m, const int n, const int k,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                           const int bs1, const ScalarType beta,
                                           /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1);
};

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
    OpA opA, OpB opB, const int m, const int n, const int k, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
    const int as0, const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
    const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1) {
  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)

  const ScalarType one(1.0), zero(0.0);

  if (beta == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, n, zero, C, cs0, cs1);
  else if (beta != one)
    KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, beta, C, cs0, cs1);

  if (alpha != zero) {
    ValueType *KOKKOS_RESTRICT pC = C;
    for (int p = 0; p < k; ++p) {
      const ValueType *KOKKOS_RESTRICT pA = A + p * as1, *KOKKOS_RESTRICT pB = B + p * bs0;
      for (int i = 0; i < m; ++i) {
        const ValueType tA(alpha * opA(pA[i * as0]));
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
        for (int j = 0; j < n; ++j) pC[i * cs0 + j * cs1] += tA * opB(pB[j * bs1]);
      }
    }
  }
  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
    OpA opA, OpB opB, const int m, const int n, const int k, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
    const int as0, const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
    const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1) {
  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)

  constexpr int mbAlgo = Algo::Gemm::Blocked::mb();
  constexpr int nbAlgo = Algo::Gemm::Blocked::mb();

  const ScalarType one(1.0), zero(0.0);

  if (beta == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, n, zero, C, cs0, cs1);
  else if (beta != one)
    KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, beta, C, cs0, cs1);

  if (alpha != zero) {
    if (m <= 0 || n <= 0 || k <= 0) return 0;
    const ValueType alpha_value(alpha);

    InnerGemmFixC<mbAlgo, nbAlgo> inner(as0, as1, bs0, bs1, cs0, cs1);
    auto gemm = [&](const int ib, const int jb, const int pb, const ValueType *KOKKOS_RESTRICT AA,
                    const ValueType *KOKKOS_RESTRICT BB,
                    /**/ ValueType *KOKKOS_RESTRICT CC) {
      const int mb = mbAlgo, nb = nbAlgo;
      for (int i = 0; i < ib; i += mb)
        for (int j = 0; j < jb; j += nb)
          inner.serial_invoke(opA, opB, alpha_value, AA + i * as0, BB + j * bs1, (i + mb) > ib ? (ib - i) : mb,
                              (j + nb) > jb ? (jb - j) : nb, pb, CC + i * cs0 + j * cs1);
    };

    const bool is_small = true;  //(m*n*k <= 64*64*64);
    if (is_small) {
      gemm(m, n, k, A, B, C);
    } else {
      // // cache blocking
      // const int
      //   nc = nb*10, kc = mb*4, mc = mb*4;

      // for (int jj=0;jj<n;jj+=nc) {
      //   const int tj = n-jj, jb = (tj < nc ? tj : nc);
      //   for (int pp=0;pp<k;pp+=kc) {
      //     const int tp = k-pp, pb = (tp < kc ? tp : kc);
      //     //const int pb = k, pp = 0;
      //     for (int ii=0;ii<m;ii+=mc) {
      //       const int ti = m-ii, ib = (ti < mc ? ti : mc);

      //       const ValueType *KOKKOS_RESTRICT AA = A+ii*as0+pp*as1;
      //       const ValueType *KOKKOS_RESTRICT BB = B+pp*bs0+jj*bs1;
      //       /**/  ValueType *KOKKOS_RESTRICT CC = C+ii*cs0+jj*cs1;

      //       gemm(ib, jb, pb, AA, BB, CC);
      //     } // for ii
      //   } // for pp
      // } // for jj
    }
  }
  return 0;
}

}  // namespace Impl

template <typename ArgAlgo>
struct [[deprecated("Use KokkosBatched::SerialGemm instead")]] SerialGemmInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const int k, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                           const ScalarType beta,
                                           /**/ ValueType *KOKKOS_RESTRICT C, const int cs0, const int cs1) {
    return Impl::SerialGemmInternal<ArgAlgo>::invoke(KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), m, n, k, alpha,
                                                     A, as0, as1, B, bs0, bs1, beta, C, cs0, cs1);
  }
};

}  // namespace KokkosBatched

#endif
