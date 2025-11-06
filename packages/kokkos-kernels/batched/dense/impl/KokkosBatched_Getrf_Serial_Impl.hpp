// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GETRF_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_GETRF_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Getrf_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AViewType, typename PivViewType>
KOKKOS_INLINE_FUNCTION static int checkGetrfInput([[maybe_unused]] const AViewType &A,
                                                  [[maybe_unused]] const PivViewType &ipiv) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::getrf: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<PivViewType>, "KokkosBatched::getrf: PivViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::getrf: AViewType must have rank 2.");
  static_assert(PivViewType::rank == 1, "KokkosBatched::getrf: PivViewType must have rank 1.");
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int m = A.extent(0), n = A.extent(1);
  const int npiv = ipiv.extent(0);
  if (npiv != Kokkos::min(m, n)) {
    Kokkos::printf(
        "KokkosBatched::getrf: the dimension of the ipiv array must "
        "satisfy ipiv.extent(0) == max(m, n): ipiv: %d, A: "
        "%d "
        "x %d \n",
        npiv, m, n);
    return 1;
  }

#endif
  return 0;
}
}  // namespace Impl

template <>
struct SerialGetrf<Algo::Getrf::Unblocked> {
  template <typename AViewType, typename PivViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &ipiv) {
    // Quick return if possible
    if (A.extent(0) == 0 || A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkGetrfInput(A, ipiv);
    if (info) return info;
    KOKKOS_IF_ON_HOST((return KokkosBatched::Impl::SerialGetrfInternalHost<Algo::Getrf::Unblocked>::invoke(A, ipiv);))
    KOKKOS_IF_ON_DEVICE(
        (return KokkosBatched::Impl::SerialGetrfInternalDevice<Algo::Getrf::Unblocked>::invoke(A, ipiv);))
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GETRF_SERIAL_IMPL_HPP_
