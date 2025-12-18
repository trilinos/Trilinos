// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PBTRS_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_PBTRS_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Pbtrs_Serial_Internal.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
namespace Impl {

template <typename AViewType, typename XViewType>
KOKKOS_INLINE_FUNCTION static int checkPbtrsInput([[maybe_unused]] const AViewType &A,
                                                  [[maybe_unused]] const XViewType &x) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::pbtrs: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::pbtrs: XViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::pbtrs: AViewType must have rank 2.");
  static_assert(XViewType::rank == 1, "KokkosBatched::pbtrs: XViewType must have rank 1.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int ldb = x.extent(0);
  const int lda = A.extent(0), n = A.extent(1);
  const int kd = lda - 1;
  if (kd < 0) {
    Kokkos::printf(
        "KokkosBatched::pbtrs: leading dimension of A must not be less than 1: %d, A: "
        "%d "
        "x %d \n",
        lda, n);
    return 1;
  }
  if (ldb < Kokkos::max(1, n)) {
    Kokkos::printf(
        "KokkosBatched::pbtrs: Dimensions of x and A do not match: x: %d, A: "
        "%d "
        "x %d \n"
        "x.extent(0) must be larger or equal to A.extent(1) \n",
        ldb, lda, n);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

//// Lower ////
template <>
struct SerialPbtrs<Uplo::Lower, Algo::Pbtrs::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x) {
    // Quick return if possible
    if (A.extent(1) == 0) return 0;
    auto info = KokkosBatched::Impl::checkPbtrsInput(A, x);
    if (info) return info;

    const int kd = A.extent(0) - 1;
    return KokkosBatched::Impl::SerialPbtrsInternalLower<Algo::Pbtrs::Unblocked>::invoke(
        A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(), x.stride(0), kd);
  }
};

//// Upper ////
template <>
struct SerialPbtrs<Uplo::Upper, Algo::Pbtrs::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x) {
    // Quick return if possible
    if (A.extent(1) == 0) return 0;
    auto info = KokkosBatched::Impl::checkPbtrsInput(A, x);
    if (info) return info;

    const int kd = A.extent(0) - 1;
    return KokkosBatched::Impl::SerialPbtrsInternalUpper<Algo::Pbtrs::Unblocked>::invoke(
        A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(), x.stride(0), kd);
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PBTRS_SERIAL_IMPL_HPP_
