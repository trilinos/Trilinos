// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_IAMAX_IMPL_HPP_
#define KOKKOSBATCHED_IAMAX_IMPL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Iamax_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========
template <typename XViewType>
KOKKOS_INLINE_FUNCTION typename XViewType::size_type SerialIamax::invoke(const XViewType &x) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::iamax: XViewType is not a Kokkos::View.");
  using size_type = typename XViewType::size_type;
  return Impl::SerialIamaxInternal::invoke(static_cast<size_type>(x.extent(0)), x.data(),
                                           static_cast<size_type>(x.stride(0)));
}

///
/// Team Impl
/// ===========
template <typename MemberType>
template <typename XViewType>
KOKKOS_INLINE_FUNCTION typename XViewType::size_type TeamIamax<MemberType>::invoke(const MemberType &member,
                                                                                   const XViewType &x) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::iamax: XViewType is not a Kokkos::View.");
  using size_type = typename XViewType::size_type;
  return Impl::TeamIamaxInternal<MemberType>::invoke(member, static_cast<size_type>(x.extent(0)), x.data(),
                                                     static_cast<size_type>(x.stride(0)));
}

///
/// TeamVector Impl
/// ===============
template <typename MemberType>
template <typename XViewType>
KOKKOS_INLINE_FUNCTION typename XViewType::size_type TeamVectorIamax<MemberType>::invoke(const MemberType &member,
                                                                                         const XViewType &x) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::iamax: XViewType is not a Kokkos::View.");
  using size_type = typename XViewType::size_type;
  return Impl::TeamVectorIamaxInternal<MemberType>::invoke(member, static_cast<size_type>(x.extent(0)), x.data(),
                                                           static_cast<size_type>(x.stride(0)));
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_IAMAX_IMPL_HPP_
