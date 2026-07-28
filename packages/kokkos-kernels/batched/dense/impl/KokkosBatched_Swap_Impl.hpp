// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SWAP_IMPL_HPP_
#define KOKKOSBATCHED_SWAP_IMPL_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Swap_Internal.hpp"

namespace KokkosBatched {
namespace Impl {

template <typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION static int checkSwapInput([[maybe_unused]] const XViewType &x,
                                                 [[maybe_unused]] const YViewType &y) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::swap: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::swap: YViewType is not a Kokkos::View.");
  static_assert(XViewType::rank() == 1 || XViewType::rank() == 2,
                "KokkosBatched::swap: XViewType must have rank 1 or 2.");
  static_assert(YViewType::rank() == 1 || YViewType::rank() == 2,
                "KokkosBatched::swap: YViewType must have rank 1 or 2.");
  static_assert(XViewType::rank() == YViewType::rank(),
                "KokkosBatched::swap: XViewType and YViewType must have the same rank.");
  static_assert(std::is_same_v<typename XViewType::value_type, typename XViewType::non_const_value_type>,
                "KokkosBatched::swap: XViewType must have non-const value type.");
  static_assert(std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>,
                "KokkosBatched::swap: YViewType must have non-const value type.");
  static_assert(swappable_elements<typename XViewType::non_const_value_type, typename YViewType::non_const_value_type>,
                "KokkosBatched::swap: XViewType and YViewType must have swappable value types.");

#ifndef NDEBUG
  if constexpr (XViewType::rank() == 1) {
    if (x.extent_int(0) != y.extent_int(0)) {
      Kokkos::printf(
          "KokkosBatched::swap: x and y must have the same length: x length "
          "= "
          "%d, y length = %d\n",
          x.extent_int(0), y.extent_int(0));
      return 1;
    }
  } else {
    const int m = x.extent_int(0);
    const int n = x.extent_int(1);
    if (m != y.extent_int(0) || n != y.extent_int(1)) {
      Kokkos::printf(
          "KokkosBatched::swap: x and y must have the same shape: x shape "
          "= (%d, %d), y shape = (%d, %d)\n",
          m, n, y.extent_int(0), y.extent_int(1));
      return 1;
    }
  }
#endif
  return 0;
}
}  // namespace Impl

///
/// Serial Impl
/// ===========
template <typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION int SerialSwap::invoke(const XViewType &x, const YViewType &y) {
  auto info = Impl::checkSwapInput(x, y);
  if (info) return info;

  if (x.size() == 0 || y.size() == 0) return 0;

  if constexpr (XViewType::rank() == 1) {
    Impl::SerialSwapInternal::invoke(x.extent_int(0), x.data(), x.stride(0), y.data(), y.stride(0));
  } else {
    Impl::SerialSwapInternal::invoke(x.extent_int(0), x.extent_int(1), x.data(), x.stride(0), x.stride(1), y.data(),
                                     y.stride(0), y.stride(1));
  }
  return 0;
}

///
/// Team Impl
/// =========

template <typename MemberType>
template <typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION int TeamSwap<MemberType>::invoke(const MemberType &member, const XViewType &x,
                                                        const YViewType &y) {
  auto info = Impl::checkSwapInput(x, y);
  if (info) return info;

  if (x.size() == 0 || y.size() == 0) return 0;

  if constexpr (XViewType::rank() == 1) {
    Impl::TeamSwapInternal<MemberType>::invoke(member, x.extent_int(0), x.data(), x.stride(0), y.data(), y.stride(0));
  } else {
    Impl::TeamSwapInternal<MemberType>::invoke(member, x.extent_int(0), x.extent_int(1), x.data(), x.stride(0),
                                               x.stride(1), y.data(), y.stride(0), y.stride(1));
  }
  return 0;
}

///
/// TeamVector Impl
/// ===============
template <typename MemberType>
template <typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorSwap<MemberType>::invoke(const MemberType &member, const XViewType &x,
                                                              const YViewType &y) {
  auto info = Impl::checkSwapInput(x, y);
  if (info) return info;

  if (x.size() == 0 || y.size() == 0) return 0;

  if constexpr (XViewType::rank() == 1) {
    Impl::TeamVectorSwapInternal<MemberType>::invoke(member, x.extent_int(0), x.data(), x.stride(0), y.data(),
                                                     y.stride(0));
  } else {
    Impl::TeamVectorSwapInternal<MemberType>::invoke(member, x.extent_int(0), x.extent_int(1), x.data(), x.stride(0),
                                                     x.stride(1), y.data(), y.stride(0), y.stride(1));
  }
  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SWAP_IMPL_HPP_
