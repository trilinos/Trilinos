// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PTTRS_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_PTTRS_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
namespace Impl {
template <typename ArgUplo, typename AlgoType>
struct SerialPttrsInternal {
  template <typename RealType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, const RealType *KOKKOS_RESTRICT d, const int ds0,
                                           const ValueType *KOKKOS_RESTRICT e, const int es0,
                                           ValueType *KOKKOS_RESTRICT b, const int bs0);
};

template <typename ArgUplo, typename AlgoType>
template <typename RealType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPttrsInternal<ArgUplo, AlgoType>::invoke(
    const int n, const RealType *KOKKOS_RESTRICT d, const int ds0, const ValueType *KOKKOS_RESTRICT e, const int es0,
    ValueType *KOKKOS_RESTRICT b, const int bs0) {
  using MayBeOpConj = std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, KokkosBlas::Impl::OpConj,
                                         KokkosBlas::Impl::OpID>;
  using OpUpper     = std::conditional_t<std::is_same_v<ArgUplo, Uplo::Upper>, MayBeOpConj, KokkosBlas::Impl::OpID>;
  using OpLower     = std::conditional_t<std::is_same_v<ArgUplo, Uplo::Lower>, MayBeOpConj, KokkosBlas::Impl::OpID>;

  OpUpper op_upper;
  OpLower op_lower;

  // Solve A * X = B using the factorization L * D * L**T or L * D * L**H
  for (int i = 1; i < n; i++) {
    b[i * bs0] -= op_upper(e[(i - 1) * es0]) * b[(i - 1) * bs0];
  }

  b[(n - 1) * bs0] /= d[(n - 1) * ds0];

  for (int i = n - 2; i >= 0; i--) {
    b[i * bs0] = b[i * bs0] / d[i * ds0] - b[(i + 1) * bs0] * op_lower(e[i * es0]);
  }

  return 0;
}
}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PTTRS_SERIAL_INTERNAL_HPP_
