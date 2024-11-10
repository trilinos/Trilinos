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
#ifndef __KOKKOSBATCHED_APPLY_Q_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_APPLY_Q_SERIAL_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_ApplyQ_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

template <>
template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
KOKKOS_INLINE_FUNCTION int SerialApplyQ<Side::Left, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(
    const AViewType &A, const tViewType &t, const BViewType &B, const wViewType &w) {
  return SerialApplyQ_LeftForwardInternal::invoke(B.extent(0), B.extent(1), A.extent(1), A.data(), A.stride_0(),
                                                  A.stride_1(), t.data(), t.stride_0(), B.data(), B.stride_0(),
                                                  B.stride_1(), w.data());
}

template <>
template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
KOKKOS_INLINE_FUNCTION int SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(
    const AViewType &A, const tViewType &t, const BViewType &B, const wViewType &w) {
  return SerialApplyQ_LeftBackwardInternal::invoke(B.extent(0), B.extent(1), A.extent(1), A.data(), A.stride_0(),
                                                   A.stride_1(), t.data(), t.stride_0(), B.data(), B.stride_0(),
                                                   B.stride_1(), w.data());
}

template <>
template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
KOKKOS_INLINE_FUNCTION int SerialApplyQ<Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(
    const AViewType &A, const tViewType &t, const BViewType &B, const wViewType &w) {
  return SerialApplyQ_RightForwardInternal::invoke(B.extent(0), B.extent(1), A.extent(1), A.data(), A.stride_0(),
                                                   A.stride_1(), t.data(), t.stride_0(), B.data(), B.stride_0(),
                                                   B.stride_1(), w.data());
}

}  // namespace KokkosBatched

#endif
