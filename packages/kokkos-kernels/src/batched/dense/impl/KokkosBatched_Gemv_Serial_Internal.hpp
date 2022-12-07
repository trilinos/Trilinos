#ifndef __KOKKOSBATCHED_GEMV_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMV_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_serial_scal_impl.hpp"
#include "KokkosBatched_InnerMultipleDotProduct_Serial_Impl.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

template <typename ArgAlgo>
struct SerialGemvInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const int m, const int n, const ScalarType alpha,
      const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
      const ValueType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
      /**/ ValueType *KOKKOS_RESTRICT y, const int ys0);
};

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialGemvInternal<Algo::Gemv::Unblocked>::invoke(
    const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    const ValueType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT y, const int ys0) {
  const ScalarType one(1.0), zero(0.0);

  // y = beta y + alpha A x
  // y (m), A(m x n), B(n)

  if (beta == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, zero, y, ys0);
  else if (beta != one)
    KokkosBlas::Impl::SerialScaleInternal::invoke(m, beta, y, ys0);

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    for (int i = 0; i < m; ++i) {
      ValueType t(0);
      const ValueType *KOKKOS_RESTRICT tA = (A + i * as0);

#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int j = 0; j < n; ++j) t += tA[j * as1] * x[j * xs0];
      y[i * ys0] += alpha * t;
    }
  }
  return 0;
}

template <>
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialGemvInternal<Algo::Gemv::Blocked>::invoke(
    const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    const ValueType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT y, const int ys0) {
  const ScalarType one(1.0), zero(0.0);

  // y = beta y + alpha A x
  // y (m), A(m x n), B(n)

  constexpr int mbAlgo = Algo::Gemv::Blocked::mb();

  if (beta == zero)
    KokkosBlas::Impl::SerialSetInternal::invoke(m, zero, y, ys0);
  else if (beta != one)
    KokkosBlas::Impl::SerialScaleInternal::invoke(m, beta, y, ys0);

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    InnerMultipleDotProduct<mbAlgo> inner(as0, as1, xs0, ys0);
    const int mb = mbAlgo;
    for (int i = 0; i < m; i += mb)
      inner.serial_invoke(alpha, A + i * as0, x, (i + mb) > m ? (m - i) : mb, n,
                          y + i * ys0);
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
