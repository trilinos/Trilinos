// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_ROT_IMPL_HPP_
#define KOKKOSBATCHED_ROT_IMPL_HPP_

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Rot_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <bool Conj, typename CType, typename SType, typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION static int checkRotInput([[maybe_unused]] const XViewType &x,
                                                [[maybe_unused]] const YViewType &y) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::rot: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::rot: YViewType is not a Kokkos::View.");
  static_assert(XViewType::rank() == 1, "KokkosBatched::rot: XViewType must have rank 1.");
  static_assert(YViewType::rank() == 1, "KokkosBatched::rot: YViewType must have rank 1.");
  static_assert(std::is_same_v<typename XViewType::value_type, typename XViewType::non_const_value_type>,
                "KokkosBatched::rot: XViewType must have non-const value type.");
  static_assert(std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>,
                "KokkosBatched::rot: YViewType must have non-const value type.");
  static_assert(!KokkosKernels::ArithTraits<CType>::is_complex, "KokkosBatched::rot: CType must be real.");
  using x_value_type = typename XViewType::non_const_value_type;
  using y_value_type = typename YViewType::non_const_value_type;
  static_assert(
      (KokkosKernels::ArithTraits<x_value_type>::is_complex && KokkosKernels::ArithTraits<y_value_type>::is_complex) ||
          (!KokkosKernels::ArithTraits<x_value_type>::is_complex &&
           !KokkosKernels::ArithTraits<y_value_type>::is_complex),
      "KokkosBatched::rot: XViewType and YViewType must be either both complex or both real.");

  if constexpr (KokkosKernels::ArithTraits<x_value_type>::is_complex) {
    if constexpr (Conj) {
      // {c,z}rot, S must be complex
      static_assert(KokkosKernels::ArithTraits<SType>::is_complex,
                    "KokkosBatched::rot: SType must be complex for complex input with Conj = true.");
    } else {
      // {cs,zd}rot, S must be real
      static_assert(!KokkosKernels::ArithTraits<SType>::is_complex,
                    "KokkosBatched::rot: SType must be real for complex input with Conj = false.");
    }
  } else {
    // {s,d} rot, S must be real
    static_assert(!KokkosKernels::ArithTraits<SType>::is_complex,
                  "KokkosBatched::rot: SType must be real for real input.");
  }

#ifndef NDEBUG
  const int n = x.extent_int(0);

  if (y.extent_int(0) != n) {
    Kokkos::printf(
        "KokkosBatched::rot: x and y must have the same length: x length "
        "= "
        "%d, y length = %d\n",
        n, y.extent_int(0));
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

///
/// Serial Impl
/// ===========

// {s,d,cs,zd}rot interface for Conj = false
// x(i) := c*x(i) + s*y(i)
// y(i) := c*y(i) - s*x(i)
//
// {c,z}rot interface for Conj = true
// x(i) := c*x(i) + s*y(i)
// y(i) := c*y(i) - conj(s)*x(i)
template <bool Conj>
template <typename XViewType, typename YViewType, typename CType, typename SType>
KOKKOS_INLINE_FUNCTION int SerialRot<Conj>::invoke(const XViewType &x, const YViewType &y, const CType c,
                                                   const SType s) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0) return 0;

  auto info = Impl::checkRotInput<Conj, CType, SType>(x, y);
  if (info) return info;

  using op = std::conditional_t<Conj, KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;
  return Impl::SerialRotInternal::invoke(op(), n, x.data(), x.stride(0), y.data(), y.stride(0), c, s);
}

///
/// Team Impl
/// ===========

// {s,d,cs,zd}rot interface for Conj = false
// x(i) := c*x(i) + s*y(i)
// y(i) := c*y(i) - s*x(i)
//
// {c,z}rot interface for Conj = true
// x(i) := c*x(i) + s*y(i)
// y(i) := c*y(i) - conj(s)*x(i)
template <typename MemberType, bool Conj>
template <typename XViewType, typename YViewType, typename CType, typename SType>
KOKKOS_INLINE_FUNCTION int TeamRot<MemberType, Conj>::invoke(const MemberType &member, const XViewType &x,
                                                             const YViewType &y, const CType c, const SType s) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0) return 0;

  auto info = Impl::checkRotInput<Conj, CType, SType>(x, y);
  if (info) return info;

  using op = std::conditional_t<Conj, KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;

  return Impl::TeamRotInternal::invoke(member, op(), n, x.data(), x.stride(0), y.data(), y.stride(0), c, s);
}

///
/// TeamVector Impl
/// ===============

// {s,d,cs,zd}rot interface for Conj = false
// x(i) := c*x(i) + s*y(i)
// y(i) := c*y(i) - s*x(i)
//
// {c,z}rot interface for Conj = true
// x(i) := c*x(i) + s*y(i)
// y(i) := c*y(i) - conj(s)*x(i)
template <typename MemberType, bool Conj>
template <typename XViewType, typename YViewType, typename CType, typename SType>
KOKKOS_INLINE_FUNCTION int TeamVectorRot<MemberType, Conj>::invoke(const MemberType &member, const XViewType &x,
                                                                   const YViewType &y, const CType c, const SType s) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0) return 0;

  auto info = Impl::checkRotInput<Conj, CType, SType>(x, y);
  if (info) return info;

  using op = std::conditional_t<Conj, KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;

  return Impl::TeamVectorRotInternal::invoke(member, op(), n, x.data(), x.stride(0), y.data(), y.stride(0), c, s);
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_ROT_IMPL_HPP_
