// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_NRM_IMPL_HPP_
#define KOKKOSBATCHED_NRM_IMPL_HPP_

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Nrm_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename XViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION static void checkNrmInput(const XViewType & /* x */, const NormViewType & /* norm */) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::nrm: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<NormViewType>, "KokkosBatched::nrm: NormViewType is not a Kokkos::View.");
  static_assert(XViewType::rank() == 1, "KokkosBatched::nrm: XViewType must have rank 1.");
  static_assert(NormViewType::rank() == 0, "KokkosBatched::nrm: NormViewType must have rank 0.");
  static_assert(std::is_same_v<typename NormViewType::value_type, typename NormViewType::non_const_value_type>,
                "KokkosBatched::nrm: NormViewType must have non-const value type.");
  static_assert(!KokkosKernels::ArithTraits<typename NormViewType::non_const_value_type>::is_complex,
                "KokkosBatched::nrm: NormViewType must be real.");
}
}  // namespace Impl

///
/// Serial Impl
/// ===========
template <typename NrmType>
template <typename XViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION int SerialNrm<NrmType>::invoke(const XViewType &x, const NormViewType &norm) {
  Impl::checkNrmInput(x, norm);
  const int n = x.extent_int(0);
  if (n == 0) {
    norm() = KokkosKernels::ArithTraits<typename NormViewType::non_const_value_type>::zero();
    return 0;
  }
  Impl::SerialNrmInternal<NrmType>::invoke(n, x.data(), x.stride(0), norm.data());
  return 0;
}

///
/// Team Impl
/// =========

template <typename MemberType, typename NrmType>
template <typename XViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION int TeamNrm<MemberType, NrmType>::invoke(const MemberType &member, const XViewType &x,
                                                                const NormViewType &norm) {
  Impl::checkNrmInput(x, norm);
  const int n = x.extent_int(0);
  if (n == 0) {
    norm() = KokkosKernels::ArithTraits<typename NormViewType::non_const_value_type>::zero();
    return 0;
  }
  Impl::TeamNrmInternal<MemberType, NrmType>::invoke(member, n, x.data(), x.stride(0), norm.data());
  return 0;
}

///
/// TeamVector Impl
/// ===============

template <typename MemberType, typename NrmType>
template <typename XViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorNrm<MemberType, NrmType>::invoke(const MemberType &member, const XViewType &x,
                                                                      const NormViewType &norm) {
  Impl::checkNrmInput(x, norm);
  const int n = x.extent_int(0);
  if (n == 0) {
    norm() = KokkosKernels::ArithTraits<typename NormViewType::non_const_value_type>::zero();
    return 0;
  }
  Impl::TeamVectorNrmInternal<MemberType, NrmType>::invoke(member, n, x.data(), x.stride(0), norm.data());
  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_NRM_IMPL_HPP_
