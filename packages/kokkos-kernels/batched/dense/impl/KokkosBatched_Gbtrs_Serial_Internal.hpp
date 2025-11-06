// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GBTRS_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_GBTRS_SERIAL_INTERNAL_HPP_

#include <Kokkos_Swap.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBlas_util.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBatched_Tbsv.hpp>
#include <KokkosBatched_Lacgv.hpp>

namespace KokkosBatched {
namespace Impl {
template <typename ArgTrans, typename ArgAlgo>
struct SerialGbtrsInternal {
  template <typename AViewType, typename PivViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &piv, const BViewType &b, const int kl,
                                           const int ku);
};

//// Non-transpose ////
template <>
template <typename AViewType, typename PivViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialGbtrsInternal<Trans::NoTranspose, Algo::Gbtrs::Unblocked>::invoke(
    const AViewType &A, const PivViewType &piv, const BViewType &b, const int kl, const int ku) {
  const int n  = A.extent(1);
  bool lnoti   = kl > 0;
  const int kd = ku + kl + 1;
  if (lnoti) {
    for (int j = 0; j < n - 1; ++j) {
      const int lm = Kokkos::min(kl, n - j - 1);
      auto l       = piv(j);
      // If pivot index is not j, swap rows l and j in b
      if (l != j) {
        Kokkos::kokkos_swap(b(l), b(j));
      }

      // Perform a rank-1 update of the remaining part of the current column
      // (ger)
      for (int i = 0; i < lm; ++i) {
        b(j + 1 + i) = b(j + 1 + i) - A(kd + i, j) * b(j);
      }
    }
  }

  // Solve U*X = b for each right hand side, overwriting B with X.
  KokkosBatched::SerialTbsv<Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, Algo::Trsv::Unblocked>::invoke(A, b,
                                                                                                           kl + ku);

  return 0;
}

//// Transpose ////
template <>
template <typename AViewType, typename PivViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialGbtrsInternal<Trans::Transpose, Algo::Gbtrs::Unblocked>::invoke(
    const AViewType &A, const PivViewType &piv, const BViewType &b, const int kl, const int ku) {
  const int n  = A.extent(1);
  bool lnoti   = kl > 0;
  const int kd = ku + kl + 1;

  // Solve U*X = b for each right hand side, overwriting B with X.
  KokkosBatched::SerialTbsv<Uplo::Upper, Trans::Transpose, Diag::NonUnit, Algo::Tbsv::Unblocked>::invoke(A, b, kl + ku);

  if (lnoti) {
    for (int j = n - 2; j >= 0; --j) {
      const int lm = Kokkos::min(kl, n - j - 1);

      // Gemv transposed
      auto a = Kokkos::subview(b, Kokkos::pair(j + 1, j + 1 + lm));
      auto x = Kokkos::subview(A, Kokkos::pair(kd, kd + lm), j);
      auto y = Kokkos::subview(b, Kokkos::pair(j, j + lm));

      KokkosBlas::Impl::SerialGemvInternal<Algo::Gemv::Unblocked>::invoke(
          1, a.extent(0), -1.0, a.data(), a.stride(0), a.stride(0), x.data(), x.stride(0), 1.0, y.data(), y.stride(0));

      // If pivot index is not j, swap rows l and j in b
      auto l = piv(j);
      if (l != j) {
        Kokkos::kokkos_swap(b(l), b(j));
      }
    }
  }

  return 0;
}

//// Conj-Transpose ////
template <>
template <typename AViewType, typename PivViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialGbtrsInternal<Trans::ConjTranspose, Algo::Gbtrs::Unblocked>::invoke(
    const AViewType &A, const PivViewType &piv, const BViewType &b, const int kl, const int ku) {
  const int n  = A.extent(1);
  bool lnoti   = kl > 0;
  const int kd = ku + kl + 1;

  // Solve U*X = b for each right hand side, overwriting B with X.
  KokkosBatched::SerialTbsv<Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit, Algo::Tbsv::Unblocked>::invoke(A, b,
                                                                                                             kl + ku);
  if (lnoti) {
    for (int j = n - 2; j >= 0; --j) {
      const int lm = Kokkos::min(kl, n - j - 1);

      // Gemv transposed
      auto a   = Kokkos::subview(b, Kokkos::pair(j + 1, j + 1 + lm));
      auto x   = Kokkos::subview(A, Kokkos::pair(kd, kd + lm), j);
      auto y   = Kokkos::subview(b, Kokkos::pair(j, j + lm));
      auto b_j = Kokkos::subview(b, Kokkos::pair(j, j + 1));

      SerialLacgv::invoke(b_j);
      KokkosBlas::Impl::SerialGemvInternal<Algo::Gemv::Unblocked>::invoke(
          KokkosBlas::Impl::OpConj(), 1, a.extent(0), -1.0, a.data(), a.stride(0), a.stride(0), x.data(), x.stride(0),
          1.0, y.data(), y.stride(0));
      SerialLacgv::invoke(b_j);

      // If pivot index is not j, swap rows l and j in b
      auto l = piv(j);
      if (l != j) {
        Kokkos::kokkos_swap(b(l), b(j));
      }
    }
  }

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GBTRS_SERIAL_INTERNAL_HPP_
