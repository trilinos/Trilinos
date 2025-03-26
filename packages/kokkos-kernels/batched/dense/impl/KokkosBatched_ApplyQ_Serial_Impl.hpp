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
#ifndef KOKKOSBATCHED_APPLY_Q_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_APPLY_Q_SERIAL_IMPL_HPP

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
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::SerialApplyQ::invoke: AViewType must be a Kokkos::View");
  static_assert(AViewType::rank() == 2, "KokkosBatched::SerialApplyQ::invoke: AViewType must have rank 2");

  static_assert(Kokkos::is_view_v<tViewType>, "KokkosBatched::SerialApplyQ::invoke: tViewType must be a Kokkos::View");
  static_assert(tViewType::rank() == 1, "KokkosBatched::SerialApplyQ::invoke: tViewType must have rank 1");

  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::SerialApplyQ::invoke: BViewType must be a Kokkos::View");
  static_assert(BViewType::rank() == 2, "KokkosBatched::SerialApplyQ::invoke: BViewType must have rank 2");

  static_assert(Kokkos::is_view_v<wViewType>, "KokkosBatched::SerialApplyQ::invoke: wViewType must be a Kokkos::View");
  static_assert(wViewType::rank() == 1, "KokkosBatched::SerialApplyQ::invoke: wViewType must have rank 1");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  if (!w.span_is_contiguous()) {
    Kokkos::printf("KokkosBatched::SerialApplyQ::invoke: w must have a contiguous span.");
    return 1;
  }
  if (A.extent_int(1) != t.extent_int(0)) {
    Kokkos::printf("KokkosBatched::SerialApplyQ::invoke: A.extent(1) is different from t.extent(0).");
    return 1;
  }
#endif

  return SerialApplyQ_LeftForwardInternal::invoke(B.extent(0), B.extent(1), A.extent(1), A.data(), A.stride_0(),
                                                  A.stride_1(), t.data(), t.stride_0(), B.data(), B.stride_0(),
                                                  B.stride_1(), w.data());
}

template <>
template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
KOKKOS_INLINE_FUNCTION int SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(
    const AViewType &A, const tViewType &t, const BViewType &B, const wViewType &w) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::SerialApplyQ::invoke: AViewType must be a Kokkos::View");
  static_assert(AViewType::rank() == 2, "KokkosBatched::SerialApplyQ::invoke: AViewType must have rank 2");

  static_assert(Kokkos::is_view_v<tViewType>, "KokkosBatched::SerialApplyQ::invoke: tViewType must be a Kokkos::View");
  static_assert(tViewType::rank() == 1, "KokkosBatched::SerialApplyQ::invoke: tViewType must have rank 1");

  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::SerialApplyQ::invoke: BViewType must be a Kokkos::View");
  static_assert(BViewType::rank() == 2, "KokkosBatched::SerialApplyQ::invoke: BViewType must have rank 2");

  static_assert(Kokkos::is_view_v<wViewType>, "KokkosBatched::SerialApplyQ::invoke: wViewType must be a Kokkos::View");
  static_assert(wViewType::rank() == 1, "KokkosBatched::SerialApplyQ::invoke: wViewType must have rank 1");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  if (!w.span_is_contiguous()) {
    Kokkos::printf("KokkosBatched::SerialApplyQ::invoke: w must have a contiguous span.");
    return 1;
  }
  if (A.extent_int(1) != t.extent_int(0)) {
    Kokkos::printf("KokkosBatched::SerialApplyQ::invoke: A.extent(1) is different from t.extent(0).");
    return 1;
  }
#endif

  return SerialApplyQ_LeftBackwardInternal::invoke(B.extent(0), B.extent(1), A.extent(1), A.data(), A.stride_0(),
                                                   A.stride_1(), t.data(), t.stride_0(), B.data(), B.stride_0(),
                                                   B.stride_1(), w.data());
}

template <>
template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
KOKKOS_INLINE_FUNCTION int SerialApplyQ<Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(
    const AViewType &A, const tViewType &t, const BViewType &B, const wViewType &w) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::SerialApplyQ::invoke: AViewType must be a Kokkos::View");
  static_assert(AViewType::rank() == 2, "KokkosBatched::SerialApplyQ::invoke: AViewType must have rank 2");

  static_assert(Kokkos::is_view_v<tViewType>, "KokkosBatched::SerialApplyQ::invoke: tViewType must be a Kokkos::View");
  static_assert(tViewType::rank() == 1, "KokkosBatched::SerialApplyQ::invoke: tViewType must have rank 1");

  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::SerialApplyQ::invoke: BViewType must be a Kokkos::View");
  static_assert(BViewType::rank() == 2, "KokkosBatched::SerialApplyQ::invoke: BViewType must have rank 2");

  static_assert(Kokkos::is_view_v<wViewType>, "KokkosBatched::SerialApplyQ::invoke: wViewType must be a Kokkos::View");
  static_assert(wViewType::rank() == 1, "KokkosBatched::SerialApplyQ::invoke: wViewType must have rank 1");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  if (!w.span_is_contiguous()) {
    Kokkos::printf("KokkosBatched::SerialApplyQ::invoke: w must have a contiguous span.");
    return 1;
  }
  if (A.extent_int(1) != t.extent_int(0)) {
    Kokkos::printf("KokkosBatched::SerialApplyQ::invoke: A.extent(1) is different from t.extent(0).");
    return 1;
  }
#endif

  return SerialApplyQ_RightForwardInternal::invoke(B.extent(0), B.extent(1), A.extent(1), A.data(), A.stride_0(),
                                                   A.stride_1(), t.data(), t.stride_0(), B.data(), B.stride_0(),
                                                   B.stride_1(), w.data());
}

}  // namespace KokkosBatched

#endif
