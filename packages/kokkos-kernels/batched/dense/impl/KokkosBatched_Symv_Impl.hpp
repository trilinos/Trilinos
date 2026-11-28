// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SYMV_IMPL_HPP_
#define KOKKOSBATCHED_SYMV_IMPL_HPP_

#include <concepts>
#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Symv_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AViewType, typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION static int checkSymvInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const XViewType &x,
                                                 [[maybe_unused]] const YViewType &y) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::symv: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::symv: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::symv: YViewType is not a Kokkos::View.");

  static_assert(AViewType::rank == 2, "KokkosBatched::symv: AViewType must have rank 2.");
  static_assert(XViewType::rank == 1, "KokkosBatched::symv: XViewType must have rank 1.");
  static_assert(YViewType::rank == 1, "KokkosBatched::symv: YViewType must have rank 1.");
#ifndef NDEBUG
  const int lda = A.extent_int(0), n = A.extent_int(1);
  const int x_len = x.extent_int(0), y_len = y.extent_int(0);
  if (x_len != y_len) {
    Kokkos::printf(
        "KokkosBatched::symv: x and y must have the same length: x_len = "
        "%d, y_len = %d\n",
        x_len, y_len);
    return 1;
  }

  if (n != x_len) {
    Kokkos::printf(
        "KokkosBatched::symv: length of x must match the number of columns of A: n = "
        "%d, x_len = %d\n",
        n, x_len);
    return 1;
  }

  if (lda < Kokkos::max(1, n)) {
    Kokkos::printf(
        "KokkosBatched::symv: leading dimension of A must not be smaller than "
        "max(1, n): "
        "lda = %d, n = %d\n",
        lda, n);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

/// Serial Batched Symv:
template <typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION int SerialSymv<ArgUplo, ArgTrans>::invoke(const ScalarType alpha, const AViewType &A,
                                                                 const XViewType &x, const ScalarType beta,
                                                                 const YViewType &y) {
  // Quick return if possible
  const int n = A.extent_int(1);
  if (n == 0 || (alpha == ScalarType(0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSymvInput(A, x, y);
  if (info) return info;

  using op     = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpReal,
                                    KokkosBlas::Impl::OpID>;

  return Impl::SerialSymvInternal::invoke<ArgUplo>(op(), sym_op(), n, alpha, A.data(), A.stride(0), A.stride(1),
                                                   x.data(), x.stride(0), beta, y.data(), y.stride(0));
}

/// Team Batched Symv:
template <typename MemberType, typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION int TeamSymv<MemberType, ArgUplo, ArgTrans>::invoke(const MemberType &member,
                                                                           const ScalarType alpha, const AViewType &A,
                                                                           const XViewType &x, const ScalarType beta,
                                                                           const YViewType &y) {
  // Quick return if possible
  const int n = A.extent_int(1);
  if (n == 0 || (alpha == ScalarType(0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSymvInput(A, x, y);
  if (info) return info;

  using op     = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpReal,
                                    KokkosBlas::Impl::OpID>;

  return Impl::TeamSymvInternal::invoke<ArgUplo>(member, op(), sym_op(), n, alpha, A.data(), A.stride(0), A.stride(1),
                                                 x.data(), x.stride(0), beta, y.data(), y.stride(0));
}

/// TeamVector Batched Symv:
template <typename MemberType, typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorSymv<MemberType, ArgUplo, ArgTrans>::invoke(const MemberType &member,
                                                                                 const ScalarType alpha,
                                                                                 const AViewType &A, const XViewType &x,
                                                                                 const ScalarType beta,
                                                                                 const YViewType &y) {
  // Quick return if possible
  const int n = A.extent_int(1);
  if (n == 0 || (alpha == ScalarType(0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSymvInput(A, x, y);
  if (info) return info;

  using op     = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpReal,
                                    KokkosBlas::Impl::OpID>;

  return Impl::TeamVectorSymvInternal::invoke<ArgUplo>(member, op(), sym_op(), n, alpha, A.data(), A.stride(0),
                                                       A.stride(1), x.data(), x.stride(0), beta, y.data(), y.stride(0));
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SYMV_IMPL_HPP_
