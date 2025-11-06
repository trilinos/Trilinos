// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PTTRS_HPP_
#define KOKKOSBATCHED_PTTRS_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
/// \brief Serial Batched Pttrs:
/// Solve Ab_l x_l = b_l for all l = 0, ..., N
///   using the factorization A = U**H * D * U or A = L * D * L**H computed by
///   Pttrf.
///
/// \tparam ArgUplo: Type specifying the form of the factorization and whether the vector E is the superdiagonal
/// of the upper bidiagonal factor U or the subdiagonal of the lower bidiagonal factor L.
/// Used only for a complex matrix A
/// The matrix has the form
///    A = U**H * D * U,  if ArgUplo = KokkosBatched::Uplo::Upper, or
///    A = L * D * L**H,  if ArgUplo = KokkosBatched::Uplo::Lower,
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Pttrs::Blocked) or unblocked
/// (KokkosBatched::Algo::Pttrs::Unblocked) algorithm to be used
///
/// \tparam DViewType: Input type for the a diagonal matrix, needs to be a 1D view
/// \tparam EViewType: Input type for the a upper/lower diagonal matrix, needs to be a 1D view
/// \tparam BViewType: Input type for the right-hand side and the solution, needs to be a 1D view
///
/// \param d [in]: n diagonal elements of the diagonal matrix D
/// \param e [in]: n-1 upper/lower diagonal elements of the diagonal matrix E
/// \param b [inout]: right-hand side and the solution, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgUplo, typename ArgAlgo>
struct SerialPttrs {
  static_assert(std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
                "KokkosBatched::pttrs: Use Uplo::Upper (or Lower) if vector E specifies the superdiagonal (or "
                "subdiagonal) of a unit bidiagonal matrix U (or L)");
  static_assert(std::is_same_v<ArgAlgo, Algo::Pttrs::Unblocked>, "KokkosBatched::pttrs: Use Algo::Pttrs::Unblocked");
  template <typename DViewType, typename EViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DViewType &d, const EViewType &e, const BViewType &b);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Pttrs_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_PTTRS_HPP_
