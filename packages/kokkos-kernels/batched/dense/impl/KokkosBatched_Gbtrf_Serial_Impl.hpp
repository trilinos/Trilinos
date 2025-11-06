// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GBTRF_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_GBTRF_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Gbtrf_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename ABViewType, typename PivViewType>
KOKKOS_INLINE_FUNCTION static int checkGbtrfInput([[maybe_unused]] const ABViewType &AB,
                                                  [[maybe_unused]] const PivViewType &ipiv,
                                                  [[maybe_unused]] const int kl, [[maybe_unused]] const int ku,
                                                  [[maybe_unused]] const int m) {
  static_assert(Kokkos::is_view_v<ABViewType>, "KokkosBatched::gbtrf: ABViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<PivViewType>, "KokkosBatched::gbtrf: PivViewType is not a Kokkos::View.");
  static_assert(ABViewType::rank == 2, "KokkosBatched::gbtrf: ABViewType must have rank 2.");
  static_assert(PivViewType::rank == 1, "KokkosBatched::gbtrf: PivViewType must have rank 1.");
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int n    = AB.extent(1);
  const int npiv = ipiv.extent(0);
  if (npiv != Kokkos::min(m, n)) {
    Kokkos::printf(
        "KokkosBatched::gbtrf: the dimension of the ipiv array must "
        "satisfy ipiv.extent(0) == max(m, n): ipiv: %d, A: "
        "%d "
        "x %d \n",
        npiv, m, n);
    return 1;
  }
  if (m < 0) {
    Kokkos::printf(
        "KokkosBatched::gbtrf: input parameter m must not be less than 0: m "
        "= "
        "%d\n",
        m);
    return 1;
  }

  if (kl < 0) {
    Kokkos::printf(
        "KokkosBatched::gbtrf: input parameter kl must not be less than 0: kl "
        "= "
        "%d\n",
        kl);
    return 1;
  }

  if (ku < 0) {
    Kokkos::printf(
        "KokkosBatched::gbtrf: input parameter ku must not be less than 0: ku "
        "= "
        "%d\n",
        ku);
    return 1;
  }

  const int lda = AB.extent(0);
  if (lda < (2 * kl + ku + 1)) {
    Kokkos::printf(
        "KokkosBatched::gbtrs: leading dimension of A must be smaller than 2 * "
        "kl + ku + 1: "
        "lda = %d, kl = %d, ku = %d\n",
        lda, kl, ku);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

template <>
struct SerialGbtrf<Algo::Gbtrf::Unblocked> {
  template <typename ABViewType, typename PivViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &AB, const PivViewType &piv, const int kl, const int ku,
                                           const int m = -1) {
    // default: m == n
    int n     = AB.extent(1);
    int m_tmp = m > 0 ? m : n;

    auto info = Impl::checkGbtrfInput(AB, piv, kl, ku, m_tmp);
    if (info) return info;

    // Quick return if possible
    if (m_tmp == 0 || n == 0) return 0;

    return Impl::SerialGbtrfInternal<Algo::Gbtrf::Unblocked>::invoke(AB, piv, kl, ku, m_tmp);
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GBTRF_SERIAL_IMPL_HPP_
