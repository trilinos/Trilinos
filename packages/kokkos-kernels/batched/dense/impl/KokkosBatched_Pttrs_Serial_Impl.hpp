// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PTTRS_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_PTTRS_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include <KokkosBlas1_scal.hpp>
#include "KokkosBatched_Pttrs_Serial_Internal.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
namespace Impl {
template <typename DViewType, typename EViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION static int checkPttrsInput([[maybe_unused]] const DViewType &d,
                                                  [[maybe_unused]] const EViewType &e,
                                                  [[maybe_unused]] const BViewType &b) {
  static_assert(Kokkos::is_view_v<DViewType>, "KokkosBatched::pttrs: DViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<EViewType>, "KokkosBatched::pttrs: EViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::pttrs: BViewType is not a Kokkos::View.");

  static_assert(DViewType::rank() == 1, "KokkosBatched::pttrs: DViewType must have rank 1.");
  static_assert(EViewType::rank() == 1, "KokkosBatched::pttrs: EViewType must have rank 1.");
  static_assert(BViewType::rank() == 1, "KokkosBatched::pttrs: BViewType must have rank 1.");

  static_assert(std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type>,
                "KokkosBatched::pttrs: BViewType must have non-const value type.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int nd  = d.extent_int(0);
  const int ne  = e.extent_int(0);
  const int ldb = b.extent_int(0);

  if ((nd > 0 || ne > 0) && (ne + 1 != nd)) {
    Kokkos::printf(
        "KokkosBatched::pttrs: Dimensions of d and e do not match: d: %d, e: %d \n"
        "e.extent(0) must be equal to d.extent(0) - 1\n",
        nd, ne);
    return 1;
  }

  if (ldb < nd) {
    Kokkos::printf(
        "KokkosBatched::pttrs: leading dimension of b must not be smaller than "
        "n: ldb = %d, n = %d\n",
        ldb, nd);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

template <typename ArgUplo>
struct SerialPttrs<ArgUplo, Algo::Pttrs::Unblocked> {
  template <typename DViewType, typename EViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DViewType &d, const EViewType &e, const BViewType &b) {
    auto info = Impl::checkPttrsInput(d, e, b);
    if (info) return info;

    // Quick return if possible
    int n = d.extent_int(0);
    if (n == 0) return 0;

    if (n == 1) {
      using ScalarType       = typename DViewType::non_const_value_type;
      const ScalarType alpha = 1.0 / d(0);
      return KokkosBlas::SerialScale::invoke(alpha, b);
    }

    // Solve A * X = B using the factorization A = L*D*L**T,
    // overwriting each right hand side vector with its solution.
    return Impl::SerialPttrsInternal<ArgUplo, Algo::Pttrs::Unblocked>::invoke(n, d.data(), d.stride(0), e.data(),
                                                                              e.stride(0), b.data(), b.stride(0));
  }
};
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PTTRS_SERIAL_IMPL_HPP_
