// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_PBTRS_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_PBTRS_SERIAL_INTERNAL_HPP_

#include "KokkosBlas_util.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Tbsv_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

///
/// Lower
///

template <typename AlgoType>
struct SerialPbtrsInternalLower {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int an, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int kd);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrsInternalLower<Algo::Pbtrs::Unblocked>::invoke(const int an,
                                                                                    const ValueType *KOKKOS_RESTRICT A,
                                                                                    const int as0, const int as1,
                                                                                    /**/ ValueType *KOKKOS_RESTRICT x,
                                                                                    const int xs0, const int kd) {
  // Solve L*X = B, overwriting B with X.
  SerialTbsvInternalLower<Algo::Tbsv::Unblocked>::invoke(false, an, A, as0, as1, x, xs0, kd);

  // Solve L**T *X = B, overwriting B with X.
  using op = std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  SerialTbsvInternalLowerTranspose<Algo::Tbsv::Unblocked>::invoke(op(), false, an, A, as0, as1, x, xs0, kd);

  return 0;
}

///
/// Upper
///

template <typename AlgoType>
struct SerialPbtrsInternalUpper {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int an, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int kd);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPbtrsInternalUpper<Algo::Pbtrs::Unblocked>::invoke(const int an,
                                                                                    const ValueType *KOKKOS_RESTRICT A,
                                                                                    const int as0, const int as1,
                                                                                    /**/ ValueType *KOKKOS_RESTRICT x,
                                                                                    const int xs0, const int kd) {
  // Solve U**T *X = B, overwriting B with X.
  using op = std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  SerialTbsvInternalUpperTranspose<Algo::Tbsv::Unblocked>::invoke(op(), false, an, A, as0, as1, x, xs0, kd);

  // Solve U*X = B, overwriting B with X.
  SerialTbsvInternalUpper<Algo::Tbsv::Unblocked>::invoke(false, an, A, as0, as1, x, xs0, kd);

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PBTRS_SERIAL_INTERNAL_HPP_
