// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GER_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_GER_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Ger_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename XViewType, typename YViewType, typename AViewType>
KOKKOS_INLINE_FUNCTION static int checkGerInput([[maybe_unused]] const XViewType &x,
                                                [[maybe_unused]] const YViewType &y,
                                                [[maybe_unused]] const AViewType &A) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::ger: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::ger: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::ger: AViewType is not a Kokkos::View.");
  static_assert(XViewType::rank == 1, "KokkosBatched::ger: XViewType must have rank 1.");
  static_assert(YViewType::rank == 1, "KokkosBatched::ger: YViewType must have rank 1.");
  static_assert(AViewType::rank == 2, "KokkosBatched::ger: AViewType must have rank 2.");
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int lda = A.extent_int(0), n = A.extent_int(1);
  const int m = x.extent_int(0);
  if (m < 0) {
    Kokkos::printf(
        "KokkosBatched::ger: input parameter m must not be less than 0: m "
        "= "
        "%d\n",
        m);
    return 1;
  }

  if (n < 0) {
    Kokkos::printf(
        "KokkosBatched::ger: input parameter n must not be less than 0: n "
        "= "
        "%d\n",
        n);
    return 1;
  }

  if (y.extent_int(0) != n) {
    Kokkos::printf(
        "KokkosBatched::ger: y must contain n elements: n "
        "= "
        "%d\n",
        n);
    return 1;
  }

  if (lda < Kokkos::max(1, m)) {
    Kokkos::printf(
        "KokkosBatched::ger: leading dimension of A must not be smaller than "
        "max(1, m): "
        "lda = %d, m = %d\n",
        lda, m);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

// T
// A: alpha * x * y**T + A
template <>
struct SerialGer<Trans::Transpose> {
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &A) {
    // Quick return if possible
    const int m = A.extent_int(0), n = A.extent_int(1);
    if (m == 0 || n == 0 || (alpha == ScalarType(0))) return 0;

    auto info = Impl::checkGerInput(x, y, A);
    if (info) return info;

    return Impl::SerialGerInternal::invoke(KokkosBlas::Impl::OpID(), m, n, alpha, x.data(), x.stride(0), y.data(),
                                           y.stride(0), A.data(), A.stride(0), A.stride(1));
  }
};

// C
// A: alpha * x * y**H + A
template <>
struct SerialGer<Trans::ConjTranspose> {
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &A) {
    // Quick return if possible
    const int m = A.extent_int(0), n = A.extent_int(1);
    if (m == 0 || n == 0 || (alpha == ScalarType(0))) return 0;

    auto info = Impl::checkGerInput(x, y, A);
    if (info) return info;

    return Impl::SerialGerInternal::invoke(KokkosBlas::Impl::OpConj(), m, n, alpha, x.data(), x.stride(0), y.data(),
                                           y.stride(0), A.data(), A.stride(0), A.stride(1));
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GER_SERIAL_IMPL_HPP_
