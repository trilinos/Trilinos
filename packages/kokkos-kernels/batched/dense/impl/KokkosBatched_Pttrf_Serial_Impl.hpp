// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PTTRF_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_PTTRF_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Pttrf_Serial_Internal.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
namespace Impl {
template <typename DViewType, typename EViewType>
KOKKOS_INLINE_FUNCTION static int checkPttrfInput([[maybe_unused]] const DViewType &d,
                                                  [[maybe_unused]] const EViewType &e) {
  static_assert(Kokkos::is_view<DViewType>::value, "KokkosBatched::pttrf: DViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<EViewType>::value, "KokkosBatched::pttrf: EViewType is not a Kokkos::View.");

  static_assert(DViewType::rank == 1, "KokkosBatched::pttrf: DViewType must have rank 1.");
  static_assert(EViewType::rank == 1, "KokkosBatched::pttrf: EViewType must have rank 1.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int nd = d.extent(0);
  const int ne = e.extent(0);

  if (ne + 1 != nd) {
    Kokkos::printf(
        "KokkosBatched::pttrf: Dimensions of d and e do not match: d: %d, e: "
        "%d \n"
        "e.extent(0) must be equal to d.extent(0) - 1\n",
        nd, ne);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

template <>
struct SerialPttrf<Algo::Pttrf::Unblocked> {
  template <typename DViewType, typename EViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DViewType &d, const EViewType &e) {
    using ScalarType = typename DViewType::non_const_value_type;
    // Quick return if possible
    if (d.extent(0) == 0) return 0;
    if (d.extent(0) == 1) return (d(0) < KokkosKernels::ArithTraits<ScalarType>::zero() ? 1 : 0);

    auto info = Impl::checkPttrfInput(d, e);
    if (info) return info;

    return Impl::SerialPttrfInternal<Algo::Pttrf::Unblocked>::invoke(d.extent(0), d.data(), d.stride(0), e.data(),
                                                                     e.stride(0));
  }
};
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PTTRF_SERIAL_IMPL_HPP_
