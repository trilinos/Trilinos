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

#ifndef KOKKOSBATCHED_PBTRS_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_PBTRS_SERIAL_INTERNAL_HPP_

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
  constexpr bool do_conj = Kokkos::ArithTraits<ValueType>::is_complex;
  SerialTbsvInternalLowerTranspose<Algo::Tbsv::Unblocked>::invoke(false, do_conj, an, A, as0, as1, x, xs0, kd);

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
  constexpr bool do_conj = Kokkos::ArithTraits<ValueType>::is_complex;
  SerialTbsvInternalUpperTranspose<Algo::Tbsv::Unblocked>::invoke(false, do_conj, an, A, as0, as1, x, xs0, kd);

  // Solve U*X = B, overwriting B with X.
  SerialTbsvInternalUpper<Algo::Tbsv::Unblocked>::invoke(false, an, A, as0, as1, x, xs0, kd);

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PBTRS_SERIAL_INTERNAL_HPP_
