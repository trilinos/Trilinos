// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_TRMV_IMPL_HPP_
#define KOKKOSBATCHED_TRMV_IMPL_HPP_

#include <concepts>
#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Trmv_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AViewType, typename XViewType>
KOKKOS_INLINE_FUNCTION static int checkTrmvInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const XViewType &x) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::trmv: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::trmv: XViewType is not a Kokkos::View.");
  static_assert(std::same_as<typename XViewType::value_type, typename XViewType::non_const_value_type>,
                "KokkosBatched::trmv: XViewType must have a non-const value type.");
  static_assert(AViewType::rank == 2, "KokkosBatched::trmv: AViewType must have rank 2.");
  static_assert(XViewType::rank == 1, "KokkosBatched::trmv: XViewType must have rank 1.");
#ifndef NDEBUG
  const int lda = A.extent_int(0), n = A.extent_int(1);

  if (n < 0) {
    Kokkos::printf(
        "KokkosBatched::trmv: input parameter n must not be less than 0: n "
        "= "
        "%d\n",
        n);
    return 1;
  }

  if (x.extent_int(0) != n) {
    Kokkos::printf(
        "KokkosBatched::trmv: x must contain n elements: n "
        "= "
        "%d\n",
        n);
    return 1;
  }

  if (lda < Kokkos::max(1, n)) {
    Kokkos::printf(
        "KokkosBatched::trmv: leading dimension of A must not be smaller than "
        "max(1, n): "
        "lda = %d, n = %d\n",
        lda, n);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

/// Batched Serial Trmv:
template <typename ArgUplo, typename ArgTrans, typename ArgDiag>
template <typename ScalarType, typename AViewType, typename XViewType>
KOKKOS_INLINE_FUNCTION int SerialTrmv<ArgUplo, ArgTrans, ArgDiag>::invoke(const ScalarType alpha, const AViewType &A,
                                                                          const XViewType &x) {
  // Quick return if possible
  const int n = A.extent_int(1);
  if (n == 0 || alpha == ScalarType(0)) return 0;

  auto info = Impl::checkTrmvInput(A, x);
  if (info) return info;

  if constexpr (std::same_as<ArgTrans, Trans::NoTranspose>) {
    return Impl::SerialTrmvInternal::invoke<ArgUplo, ArgDiag>(n, alpha, A.data(), A.stride(0), A.stride(1), x.data(),
                                                              x.stride(0));
  } else {
    using op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                  KokkosBlas::Impl::OpID>;

    return Impl::SerialTrmvInternal::invoke<ArgUplo, ArgDiag>(op(), n, alpha, A.data(), A.stride(0), A.stride(1),
                                                              x.data(), x.stride(0));
  }
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_TRMV_IMPL_HPP_
