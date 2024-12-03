//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef KOKKOSBATCHED_PTTRS_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_PTTRS_SERIAL_IMPL_HPP_

#include <KokkosBatched_Util.hpp>
#include <KokkosBlas1_scal.hpp>
#include "KokkosBatched_Pttrs_Serial_Internal.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

template <typename DViewType, typename EViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION static int checkPttrsInput([[maybe_unused]] const DViewType &d,
                                                  [[maybe_unused]] const EViewType &e,
                                                  [[maybe_unused]] const BViewType &b) {
  static_assert(Kokkos::is_view_v<DViewType>, "KokkosBatched::pttrs: DViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<EViewType>, "KokkosBatched::pttrs: EViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::pttrs: BViewType is not a Kokkos::View.");

  static_assert(DViewType::rank == 1, "KokkosBatched::pttrs: DViewType must have rank 1.");
  static_assert(EViewType::rank == 1, "KokkosBatched::pttrs: EViewType must have rank 1.");
  static_assert(BViewType::rank == 1, "KokkosBatched::pttrs: BViewType must have rank 1.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int nd  = d.extent(0);
  const int ne  = e.extent(0);
  const int ldb = b.extent(0);

  if (ne + 1 != nd) {
    Kokkos::printf(
        "KokkosBatched::pttrs: Dimensions of d and e do not match: d: %d, e: "
        "%d \n"
        "e.extent(0) must be equal to d.extent(0) - 1\n",
        nd, ne);
    return 1;
  }

  if (ldb < Kokkos::max(1, nd)) {
    Kokkos::printf(
        "KokkosBatched::pttrs: Dimensions of d and b do not match: d: %d, b: "
        "%d \n"
        "b.extent(0) must be larger or equal to d.extent(0) \n",
        ldb, nd);
    return 1;
  }
#endif
  return 0;
}

template <typename ArgUplo>
struct SerialPttrs<ArgUplo, Algo::Pttrs::Unblocked> {
  template <typename DViewType, typename EViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DViewType &d, const EViewType &e, const BViewType &b) {
    // Quick return if possible
    if (d.extent(0) == 0) return 0;

    auto info = checkPttrsInput(d, e, b);
    if (info) return info;

    using ScalarType = typename DViewType::non_const_value_type;
    int n            = d.extent(0);

    if (n == 1) {
      const ScalarType alpha = 1.0 / d(0);
      return KokkosBlas::SerialScale::invoke(alpha, b);
    }

    // Solve A * X = B using the factorization A = L*D*L**T,
    // overwriting each right hand side vector with its solution.
    return SerialPttrsInternal<ArgUplo, Algo::Pttrs::Unblocked>::invoke(n, d.data(), d.stride(0), e.data(), e.stride(0),
                                                                        b.data(), b.stride(0));
  }
};
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PTTRS_SERIAL_IMPL_HPP_
