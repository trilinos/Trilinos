// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_LASWP_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_LASWP_SERIAL_IMPL_HPP_

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Laswp_Serial_Internal.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
namespace Impl {

template <typename PivViewType, typename AViewType>
KOKKOS_INLINE_FUNCTION static int checkLaswpInput(const PivViewType &piv, const AViewType &A) {
  static_assert(Kokkos::is_view_v<PivViewType>, "KokkosBatched::laswp: PivViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::laswp: AViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 1 || AViewType::rank == 2, "KokkosBatched::laswp: AViewType must have rank 1 or 2.");
  static_assert(PivViewType::rank == 1, "KokkosBatched::laswp: PivViewType must have rank 1.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int npiv = piv.extent(0);
  const int lda  = A.extent(0);
  if (npiv > lda) {
    Kokkos::printf(
        "KokkosBatched::laswp: the dimension of the ipiv array must "
        "satisfy ipiv.extent(0) <= A.extent(0): ipiv: %d, A: "
        "%d \n",
        npiv, lda);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

///
/// Serial Internal Impl
/// ========================

///
//// Forward pivot apply
///

template <>
struct SerialLaswp<Direct::Forward> {
  template <typename PivViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const PivViewType &piv, const AViewType &A) {
    auto info = KokkosBatched::Impl::checkLaswpInput(piv, A);
    if (info) return info;

    if constexpr (AViewType::rank == 1) {
      const int plen = piv.extent(0), ps0 = piv.stride(0), as0 = A.stride(0);
      return KokkosBatched::Impl::SerialLaswpVectorForwardInternal::invoke(plen, piv.data(), ps0, A.data(), as0);
    } else if constexpr (AViewType::rank == 2) {
      // row permutation
      const int plen = piv.extent(0), ps0 = piv.stride(0), n = A.extent(1), as0 = A.stride(0), as1 = A.stride(1);
      return KokkosBatched::Impl::SerialLaswpMatrixForwardInternal::invoke(n, plen, piv.data(), ps0, A.data(), as0,
                                                                           as1);
    }
    return 0;
  }
};

///
/// Backward pivot apply
///

template <>
struct SerialLaswp<Direct::Backward> {
  template <typename PivViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const PivViewType &piv, const AViewType &A) {
    auto info = KokkosBatched::Impl::checkLaswpInput(piv, A);
    if (info) return info;

    if constexpr (AViewType::rank == 1) {
      const int plen = piv.extent(0), ps0 = piv.stride(0), as0 = A.stride(0);
      return KokkosBatched::Impl::SerialLaswpVectorBackwardInternal::invoke(plen, piv.data(), ps0, A.data(), as0);
    } else if constexpr (AViewType::rank == 2) {
      // row permutation
      const int plen = piv.extent(0), ps0 = piv.stride(0), n = A.extent(1), as0 = A.stride(0), as1 = A.stride(1);
      return KokkosBatched::Impl::SerialLaswpMatrixBackwardInternal::invoke(n, plen, piv.data(), ps0, A.data(), as0,
                                                                            as1);
    }
    return 0;
  }
};
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_LASWP_SERIAL_IMPL_HPP_
