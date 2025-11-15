// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GETRS_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_GETRS_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Getrs_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION static int checkGetrsInput([[maybe_unused]] const AViewType &A,
                                                  [[maybe_unused]] const BViewType &b) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::getrs: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::getrs: BViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::getrs: AViewType must have rank 2.");
  static_assert(BViewType::rank == 1, "KokkosBatched::getrs: BViewType must have rank 1.");
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int lda = A.extent(0), n = A.extent(1);
  if (lda < Kokkos::max(1, n)) {
    Kokkos::printf(
        "KokkosBatched::getrs: the leading dimension of the array A must "
        "satisfy lda >= max(1, n): A: "
        "%d "
        "x %d \n",
        lda, n);
    return 1;
  }

  const int ldb = b.extent(0);
  if (ldb < Kokkos::max(1, n)) {
    Kokkos::printf(
        "KokkosBatched::getrs: the leading dimension of the array b must "
        "satisfy ldb >= max(1, n): b: %d, A: "
        "%d "
        "x %d \n",
        ldb, lda, n);
    return 1;
  }
#endif
  return 0;
}

}  // namespace Impl

template <typename ArgTrans>
struct SerialGetrs<ArgTrans, Algo::Getrs::Unblocked> {
  template <typename AViewType, typename PivViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &piv, const BViewType &b) {
    // quick return if possible
    if (A.extent(1) == 0) return 0;

    auto info = Impl::checkGetrsInput(A, b);
    if (info) return info;

    return Impl::SerialGetrsInternal<ArgTrans, Algo::Getrs::Unblocked>::invoke(A, piv, b);
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GETRF_SERIAL_IMPL_HPP_
