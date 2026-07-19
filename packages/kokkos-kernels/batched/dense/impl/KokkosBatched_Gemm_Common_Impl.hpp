// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMM_COMMON_IMPL_HPP
#define KOKKOSBATCHED_GEMM_COMMON_IMPL_HPP

#include <Kokkos_DynRankView.hpp>
#include "KokkosBlas_util.hpp"
#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename ArgTransA, typename ArgTransB, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION static int checkGemmInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const BViewType &B,
                                                 [[maybe_unused]] const CViewType &C) {
  if constexpr (Kokkos::is_view_v<AViewType>) {
    static_assert(AViewType::rank <= 2, "KokkosBatched::gemm: AViewType must have rank 0, 1 or 2.");
  } else if constexpr (Kokkos::is_dyn_rank_view_v<AViewType>) {
#ifndef NDEBUG
    if (A.rank() > 2) {
      Kokkos::printf("KokkosBatched::gemm: AViewType must have rank 0, 1 or 2. Provided rank: %d\n", A.rank());
      return 1;
    }
#endif
  } else {
    static_assert(Kokkos::is_view_v<AViewType> || Kokkos::is_dyn_rank_view_v<AViewType>,
                  "KokkosBatched::gemm: AViewType must be either a Kokkos::View or a Kokkos::DynRankView.");
  }

  if constexpr (Kokkos::is_view_v<BViewType>) {
    static_assert(BViewType::rank <= 2, "KokkosBatched::gemm: BViewType must have rank 0, 1 or 2.");
  } else if constexpr (Kokkos::is_dyn_rank_view_v<BViewType>) {
#ifndef NDEBUG
    if (B.rank() > 2) {
      Kokkos::printf("KokkosBatched::gemm: BViewType must have rank 0, 1 or 2. Provided rank: %d\n", B.rank());
      return 1;
    }
#endif
  } else {
    static_assert(Kokkos::is_view_v<BViewType> || Kokkos::is_dyn_rank_view_v<BViewType>,
                  "KokkosBatched::gemm: BViewType must be either a Kokkos::View or a Kokkos::DynRankView.");
  }

  if constexpr (Kokkos::is_view_v<CViewType>) {
    static_assert(CViewType::rank <= 2, "KokkosBatched::gemm: CViewType must have rank 0, 1 or 2.");
  } else if constexpr (Kokkos::is_dyn_rank_view_v<CViewType>) {
#ifndef NDEBUG
    if (C.rank() > 2) {
      Kokkos::printf("KokkosBatched::gemm: CViewType must have rank 0, 1 or 2. Provided rank: %d\n", C.rank());
      return 1;
    }
#endif
  } else {
    static_assert(Kokkos::is_view_v<CViewType> || Kokkos::is_dyn_rank_view_v<CViewType>,
                  "KokkosBatched::gemm: CViewType must be either a Kokkos::View or a Kokkos::DynRankView.");
  }

#ifndef NDEBUG
  const int A_extent_0 = get_extent_int(A, 0);
  const int A_extent_1 = get_extent_int(A, 1);
  const int B_extent_0 = get_extent_int(B, 0);
  const int B_extent_1 = get_extent_int(B, 1);
  const int C_extent_0 = get_extent_int(C, 0);
  const int C_extent_1 = get_extent_int(C, 1);

  const int m = C_extent_0, n = C_extent_1;
  const int lda = A_extent_0;
  const int ldb = B_extent_0;

  const int ka = std::is_same_v<ArgTransA, Trans::NoTranspose> ? A_extent_1 : A_extent_0;
  const int kb = std::is_same_v<ArgTransB, Trans::NoTranspose> ? B_extent_0 : B_extent_1;

  if (ka != kb) {
    Kokkos::printf(
        "KokkosBatched::gemm: Dimensions of A and B do not match: A: %d x %d, "
        "B: %d x %d\n",
        A_extent_0, A_extent_1, B_extent_0, B_extent_1);
    return 1;
  }

  const int nrowa = std::is_same_v<ArgTransA, Trans::NoTranspose> ? m : ka;
  const int nrowb = std::is_same_v<ArgTransB, Trans::NoTranspose> ? kb : n;

  if (lda < Kokkos::max(1, nrowa)) {
    Kokkos::printf(
        "KokkosBatched::gemm: leading dimension of A must not be smaller than "
        "max(1, nrowa): "
        "lda = %d, nrowa = %d\n",
        lda, nrowa);
    return 1;
  }
  if (ldb < Kokkos::max(1, nrowb)) {
    Kokkos::printf(
        "KokkosBatched::gemm: leading dimension of B must not be smaller than "
        "max(1, nrowb): "
        "ldb = %d, nrowb = %d\n",
        ldb, nrowb);
    return 1;
  }

#endif

  return 0;
}
}  // namespace Impl
}  // namespace KokkosBatched

#endif
