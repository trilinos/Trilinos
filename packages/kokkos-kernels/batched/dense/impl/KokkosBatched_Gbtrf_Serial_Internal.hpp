// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GBTRF_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_GBTRF_SERIAL_INTERNAL_HPP_

#include <Kokkos_Swap.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBatched_Iamax.hpp>
#include "KokkosBatched_Ger_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AlgoType>
struct SerialGbtrfInternal {
  template <typename ABViewType, typename PivViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &AB, const PivViewType &piv, const int kl, const int ku,
                                           const int m);
};

template <>
template <typename ABViewType, typename PivViewType>
KOKKOS_INLINE_FUNCTION int SerialGbtrfInternal<Algo::Gbtrf::Unblocked>::invoke(const ABViewType &AB,
                                                                               const PivViewType &piv, const int kl,
                                                                               const int ku, const int m) {
  const int n = AB.extent(1);
  // Upper bandwidth of U factor
  const int kv = ku + kl;

  // Gaussian elimination with partial pivoting
  // Set fill-in elements in columns KU+2 to KV to zero
  for (int j = ku + 1; j < Kokkos::min(kv, n); ++j) {
    for (int i = kv - j + 1; i < kl; ++i) {
      AB(i, j) = 0;
    }
  }

  // JU is the index of the last column affected by the current stage of
  // the factorization
  int ju   = 0;
  int info = 0;
  for (int j = 0; j < Kokkos::min(m, n); ++j) {
    // Set fill-in elements in column J+KV to zero
    if (j + kv < n) {
      for (int i = 0; i < kl; ++i) {
        AB(i, j + kv) = 0;
      }
    }

    // Find pivot and test for singularity. KM is the number of subdiagonals
    // elements in the current column.
    int km          = Kokkos::min(kl, m - j - 1);
    auto cur_col_AB = Kokkos::subview(AB, Kokkos::pair<int, int>(kv, kv + km + 1), j);
    int jp          = SerialIamax::invoke(cur_col_AB);
    piv(j)          = jp + j;

    if (AB(kv + jp, j) == 0) {
      // If pivot is zero, set INFO to the index of the pivot unless a
      // zero pivot has already been found.
      if (info == 0) info = j + 1;
    } else {
      ju = Kokkos::max(ju, Kokkos::min(j + ku + jp, n - 1));

      // Apply the interchange to columns J to JU
      if (jp != 0) {
        for (int k = 0; k < ju - j + 1; ++k) {
          Kokkos::kokkos_swap(AB(kv + jp - k, j + k), AB(kv - k, j + k));
        }
      }

      if (km > 0) {
        // Compute multipliers
        const auto alpha = 1.0 / AB(kv, j);
        auto sub_col_AB  = Kokkos::subview(AB, Kokkos::pair<int, int>(kv + 1, kv + km + 1), j);
        KokkosBlas::SerialScale::invoke(alpha, sub_col_AB);

        // Update trailing submatrix within the band
        if (ju > j) {
          auto x = Kokkos::subview(AB, Kokkos::pair<int, int>(kv + 1, kv + km + 1), j);

          // dger or zgeru with alpha = -1.0
          const int abs0 = AB.stride(0), abs1 = AB.stride(1);
          Impl::SerialGerInternal::invoke(KokkosBlas::Impl::OpID(), km, ju - j, -1.0, &AB(kv + 1, j), abs0,
                                          &AB(kv - 1, j + 1), (abs1 - abs0), &AB(kv, j + 1), abs0, (abs1 - abs0));
        }
      }
    }
  }
  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GBTRF_SERIAL_INTERNAL_HPP_
