// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SPR_IMPL_HPP_
#define KOKKOSBATCHED_SPR_IMPL_HPP_

#include <concepts>
#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Spr_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename XViewType, typename APViewType>
KOKKOS_INLINE_FUNCTION static int checkSprInput([[maybe_unused]] const XViewType &x,
                                                [[maybe_unused]] const APViewType &ap) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::spr: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<APViewType>, "KokkosBatched::spr: APViewType is not a Kokkos::View.");
  static_assert(XViewType::rank == 1, "KokkosBatched::spr: XViewType must have rank 1.");
  static_assert(APViewType::rank == 1, "KokkosBatched::spr: APViewType must have rank 1.");
#ifndef NDEBUG
  const int n = x.extent_int(0);

  if (ap.extent_int(0) < n * (n + 1) / 2) {
    Kokkos::printf(
        "KokkosBatched::spr: size of packed A must not be smaller than n*(n+1)/2: "
        "size = %d, n = %d\n",
        ap.extent_int(0), n);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

// Serial spr interface
// A: alpha * x * x**T + A
template <typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename XViewType, typename APViewType>
KOKKOS_INLINE_FUNCTION int SerialSpr<ArgUplo, ArgTrans>::invoke(const ScalarType alpha, const XViewType &x,
                                                                const APViewType &ap) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0 || (alpha == ScalarType(0))) return 0;
  auto info = Impl::checkSprInput(x, ap);
  if (info) return info;
  using op     = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpReal,
                                    KokkosBlas::Impl::OpID>;
  return Impl::SerialSprInternal::invoke<ArgUplo>(op(), sym_op(), n, alpha, x.data(), x.stride(0), ap.data(),
                                                  ap.stride(0));
}

// Team spr interface
// A: alpha * x * x**T + A
template <typename MemberType, typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename XViewType, typename APViewType>
KOKKOS_INLINE_FUNCTION int TeamSpr<MemberType, ArgUplo, ArgTrans>::invoke(const MemberType &member,
                                                                          const ScalarType alpha, const XViewType &x,
                                                                          const APViewType &ap) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0 || (alpha == ScalarType(0))) return 0;
  auto info = Impl::checkSprInput(x, ap);
  if (info) return info;
  using op     = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpReal,
                                    KokkosBlas::Impl::OpID>;
  return Impl::TeamSprInternal::invoke<ArgUplo>(member, op(), sym_op(), n, alpha, x.data(), x.stride(0), ap.data(),
                                                ap.stride(0));
}

// TeamVector spr interface
// A: alpha * x * x**T + A
template <typename MemberType, typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename XViewType, typename APViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorSpr<MemberType, ArgUplo, ArgTrans>::invoke(const MemberType &member,
                                                                                const ScalarType alpha,
                                                                                const XViewType &x,
                                                                                const APViewType &ap) {
  // Quick return if possible
  const int n = x.extent_int(0);
  if (n == 0 || (alpha == ScalarType(0))) return 0;
  auto info = Impl::checkSprInput(x, ap);
  if (info) return info;
  using op     = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpReal,
                                    KokkosBlas::Impl::OpID>;
  return Impl::TeamVectorSprInternal::invoke<ArgUplo>(member, op(), sym_op(), n, alpha, x.data(), x.stride(0),
                                                      ap.data(), ap.stride(0));
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SPR_IMPL_HPP_
