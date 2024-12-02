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
#ifndef KOKKOSBATCHED_PBTRF_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_PBTRF_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Pbtrf_Serial_Internal.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

template <typename ABViewType>
KOKKOS_INLINE_FUNCTION static int checkPbtrfInput([[maybe_unused]] const ABViewType &Ab) {
  static_assert(Kokkos::is_view_v<ABViewType>, "KokkosBatched::pbtrf: ABViewType is not a Kokkos::View.");
  static_assert(ABViewType::rank == 2, "KokkosBatched::pbtrf: ABViewType must have rank 2.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int kd = Ab.extent(0) - 1;
  if (kd < 0) {
    Kokkos::printf(
        "KokkosBatched::pbtrf: input parameter kd must not be less than 0: kd "
        "= "
        "%d\n",
        kd);
    return 1;
  }
#endif
  return 0;
}

//// Lower ////
template <>
struct SerialPbtrf<Uplo::Lower, Algo::Pbtrf::Unblocked> {
  template <typename ABViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &Ab) {
    // Quick return if possible
    const int n = Ab.extent(1);
    if (n == 0) return 0;

    auto info = checkPbtrfInput(Ab);
    if (info) return info;

    const int kd = Ab.extent(0) - 1;
    return SerialPbtrfInternalLower<Algo::Pbtrf::Unblocked>::invoke(n, Ab.data(), Ab.stride_0(), Ab.stride_1(), kd);
  }
};

//// Upper ////
template <>
struct SerialPbtrf<Uplo::Upper, Algo::Pbtrf::Unblocked> {
  template <typename ABViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &Ab) {
    // Quick return if possible
    const int n = Ab.extent(1);
    if (n == 0) return 0;

    auto info = checkPbtrfInput(Ab);
    if (info) return info;

    const int kd = Ab.extent(0) - 1;
    return SerialPbtrfInternalUpper<Algo::Pbtrf::Unblocked>::invoke(n, Ab.data(), Ab.stride_0(), Ab.stride_1(), kd);
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PBTRF_SERIAL_IMPL_HPP_
