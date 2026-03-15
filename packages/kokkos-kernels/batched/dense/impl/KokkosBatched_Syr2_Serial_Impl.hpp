// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SYR2_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_SYR2_SERIAL_IMPL_HPP_

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Syr2_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename XViewType, typename YViewType, typename AViewType>
KOKKOS_INLINE_FUNCTION static int checkSyr2Input([[maybe_unused]] const XViewType &x,
                                                 [[maybe_unused]] const YViewType &y,
                                                 [[maybe_unused]] const AViewType &A) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::syr2: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::syr2: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::syr2: AViewType is not a Kokkos::View.");
  static_assert(XViewType::rank() == 1, "KokkosBatched::syr2: XViewType must have rank 1.");
  static_assert(YViewType::rank() == 1, "KokkosBatched::syr2: YViewType must have rank 1.");
  static_assert(AViewType::rank() == 2, "KokkosBatched::syr2: AViewType must have rank 2.");
#ifndef NDEBUG
  const int lda = A.extent_int(0), n = A.extent_int(1);

  if (n < 0) {
    Kokkos::printf(
        "KokkosBatched::syr2: input parameter n must not be less than 0: n "
        "= "
        "%d\n",
        n);
    return 1;
  }

  if (x.extent_int(0) != n) {
    Kokkos::printf(
        "KokkosBatched::syr2: x must contain n elements: n "
        "= "
        "%d\n",
        n);
    return 1;
  }

  if (y.extent_int(0) != n) {
    Kokkos::printf(
        "KokkosBatched::syr2: y must contain n elements: n "
        "= "
        "%d\n",
        n);
    return 1;
  }

  if (lda < Kokkos::max(1, n)) {
    Kokkos::printf(
        "KokkosBatched::syr2: leading dimension of A must not be smaller than "
        "max(1, n): "
        "lda = %d, n = %d\n",
        lda, n);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

// {s,d,c,z}syr2 interface
// L T
// A: alpha * x * y**T + alpha * y * x**T + A
template <>
struct SerialSyr2<Uplo::Lower, Trans::Transpose> {
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &A) {
    // Quick return if possible
    const int n = A.extent_int(1);
    if (n == 0 || (alpha == ScalarType(0))) return 0;

    auto info = Impl::checkSyr2Input(x, y, A);
    if (info) return info;

    return Impl::SerialSyr2InternalLower::invoke(KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), n, alpha, x.data(),
                                                 x.stride(0), y.data(), y.stride(0), A.data(), A.stride(0),
                                                 A.stride(1));
  }
};

// {s,d,c,z}syr2 interface
// U T
// A: alpha * x * y**T + alpha * y * x**T + A
template <>
struct SerialSyr2<Uplo::Upper, Trans::Transpose> {
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &A) {
    // Quick return if possible
    const int n = A.extent_int(1);
    if (n == 0 || (alpha == ScalarType(0))) return 0;

    auto info = Impl::checkSyr2Input(x, y, A);
    if (info) return info;

    return Impl::SerialSyr2InternalUpper::invoke(KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), n, alpha, x.data(),
                                                 x.stride(0), y.data(), y.stride(0), A.data(), A.stride(0),
                                                 A.stride(1));
  }
};

// {c,z}her2 interface
// L C
// A: alpha * x * y**H + conjg(alpha) * y * x**H + A
template <>
struct SerialSyr2<Uplo::Lower, Trans::ConjTranspose> {
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &A) {
    // Quick return if possible
    const int n = A.extent_int(1);
    if (n == 0 || (alpha == ScalarType(0))) return 0;

    auto info = Impl::checkSyr2Input(x, y, A);
    if (info) return info;

    return Impl::SerialSyr2InternalLower::invoke(KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpReal(), n, alpha,
                                                 x.data(), x.stride(0), y.data(), y.stride(0), A.data(), A.stride(0),
                                                 A.stride(1));
  }
};

// {c,z}her2 interface
// U C
// A: alpha * x * y**H + conjg(alpha) * y * x**H + A
template <>
struct SerialSyr2<Uplo::Upper, Trans::ConjTranspose> {
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &A) {
    // Quick return if possible
    const int n = A.extent_int(1);
    if (n == 0 || (alpha == ScalarType(0))) return 0;

    auto info = Impl::checkSyr2Input(x, y, A);
    if (info) return info;

    return Impl::SerialSyr2InternalUpper::invoke(KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpReal(), n, alpha,
                                                 x.data(), x.stride(0), y.data(), y.stride(0), A.data(), A.stride(0),
                                                 A.stride(1));
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SYR2_SERIAL_IMPL_HPP_
