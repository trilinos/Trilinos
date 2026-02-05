// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GBTRS_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_GBTRS_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Gbtrs_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AViewType, typename PivViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION static int checkGbtrsInput([[maybe_unused]] const AViewType &A,
                                                  [[maybe_unused]] const PivViewType &ipiv,
                                                  [[maybe_unused]] const BViewType &b, [[maybe_unused]] const int kl,
                                                  [[maybe_unused]] const int ku) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::gbtrs: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<PivViewType>, "KokkosBatched::gbtrs: PivViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::gbtrs: BViewType is not a Kokkos::View.");
  static_assert(AViewType::rank() == 2, "KokkosBatched::gbtrs: AViewType must have rank 2.");
  static_assert(PivViewType::rank() == 1, "KokkosBatched::gbtrs: PivViewType must have rank 1.");
  static_assert(BViewType::rank() == 1, "KokkosBatched::gbtrs: BViewType must have rank 1.");

  static_assert(std::is_integral_v<typename PivViewType::non_const_value_type>,
                "KokkosBatched::gbtrs: Value type of PivViewType must be an integral type.");
  static_assert(std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type>,
                "KokkosBatched::gbtrs: BViewType must have non-const value type.");
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  if (kl < 0) {
    Kokkos::printf(
        "KokkosBatched::gbtrs: input parameter kl must not be less than 0: kl "
        "= "
        "%d\n",
        kl);
    return 1;
  }

  if (ku < 0) {
    Kokkos::printf(
        "KokkosBatched::gbtrs: input parameter ku must not be less than 0: ku "
        "= "
        "%d\n",
        ku);
    return 1;
  }

  const int lda = A.extent_int(0), n = A.extent_int(1);
  if (lda < (2 * kl + ku + 1)) {
    Kokkos::printf(
        "KokkosBatched::gbtrs: leading dimension of A must not be smaller than 2 * "
        "kl + ku + 1: "
        "lda = %d, kl = %d, ku = %d\n",
        lda, kl, ku);
    return 1;
  }

  const int ldb = b.extent_int(0);
  if (ldb < n) {
    Kokkos::printf(
        "KokkosBatched::gbtrs: leading dimension of b must not be smaller than n: "
        "ldb = %d, n = %d\n",
        ldb, n);
    return 1;
  }

  const int npiv = ipiv.extent_int(0);
  if (npiv != n) {
    Kokkos::printf(
        "KokkosBatched::gbtrs: the dimension of the ipiv array must "
        "satisfy ipiv.extent(0) == n: ipiv: %d, n: %d\n",
        npiv, n);
    return 1;
  }

#endif
  return 0;
}
}  // namespace Impl

template <typename ArgTrans>
struct SerialGbtrs<ArgTrans, Algo::Gbtrs::Unblocked> {
  template <typename AViewType, typename PivViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &piv, const BViewType &b, const int kl,
                                           const int ku) {
    auto info = Impl::checkGbtrsInput(A, piv, b, kl, ku);
    if (info) return info;

    // quick return if possible
    if (A.extent(1) == 0) return 0;

    return Impl::SerialGbtrsInternal<ArgTrans, Algo::Gbtrs::Unblocked>::invoke(A, piv, b, kl, ku);
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GBTRS_SERIAL_IMPL_HPP_
