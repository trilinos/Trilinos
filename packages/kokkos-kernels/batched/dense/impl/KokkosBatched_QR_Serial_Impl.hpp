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
#ifndef KOKKOSBATCHED_QR_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_QR_SERIAL_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_QR_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

template <>
template <typename AViewType, typename tViewType, typename wViewType>
KOKKOS_INLINE_FUNCTION int SerialQR<Algo::QR::Unblocked>::invoke(const AViewType &A, const tViewType &t,
                                                                 const wViewType &w) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::SerialQR::invoke: AViewType must be a Kokkos::View");
  static_assert(AViewType::rank() == 2, "KokkosBatched::SerialQR::invoke: AViewType must have rank 2");

  static_assert(Kokkos::is_view_v<tViewType>, "KokkosBatched::SerialQR::invoke: tViewType must be a Kokkos::View");
  static_assert(tViewType::rank() == 1, "KokkosBatched::SerialQR::invoke: tViewType must have rank 1");

  static_assert(Kokkos::is_view_v<wViewType>, "KokkosBatched::SerialQR::invoke: wViewType must be a Kokkos::View");
  static_assert(wViewType::rank() == 1, "KokkosBatched::SerialQR::invoke: wViewType must have rank 1");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  if (!w.span_is_contiguous()) {
    Kokkos::printf("KokkosBatched::SerialQR::invoke: w must have a contiguous span.");
    return 1;
  }
  if (A.extent_int(1) != t.extent_int(0)) {
    Kokkos::printf("KokkosBatched::SerialQR::invoke: A.extent(1) is different from t.extent(0).");
    return 1;
  }
#endif

  return Impl::SerialQR_Internal::invoke(A.extent(0), A.extent(1), A.data(), A.stride_0(), A.stride_1(), t.data(),
                                         t.stride_0(), w.data());
}

}  // namespace KokkosBatched

#endif
