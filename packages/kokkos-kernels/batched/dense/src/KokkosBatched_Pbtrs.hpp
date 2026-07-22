// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PBTRS_HPP_
#define KOKKOSBATCHED_PBTRS_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Pbtrs:
/// Solve Ab_l x_l = b_l for all l = 0, ..., N
///   using the Cholesky factorization A = U**H * U or A = L * L**H computed by
///   Pbtrf.
/// The matrix has the form
///    A = U**H * U ,  if ArgUplo = KokkosBatched::Uplo::Upper, or
///    A = L  * L**H,  if ArgUplo = KokkosBatched::Uplo::Lower,
/// where U is an upper triangular matrix, U**H is the transpose of U, and
/// L is lower triangular matrix, L**H is the transpose of L.
///
/// \tparam ArgUplo: Type indicating whether A is the upper (Uplo::Upper) or lower (Uplo::Lower) triangular matrix
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Pbtrs::Blocked) or unblocked
/// (KokkosBatched::Algo::Pbtrs::Unblocked) algorithm to be used
///
/// \tparam ABViewType: Input type for a banded matrix, needs to be a 2D view
/// \tparam BViewType: Input type for a right-hand side and the solution, needs to be a 1D view
///
/// \param ab [in]: ab is a ldab by n banded matrix, with ( kd + 1 ) diagonals
/// \param b  [inout]: right-hand side and the solution, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgUplo, typename ArgAlgo>
struct SerialPbtrs {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::pbtrs: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::is_same_v<ArgAlgo, Algo::Pbtrs::Unblocked>, "KokkosBatched::pbtrs: Use Algo::Pbtrs::Unblocked");
  template <typename ABViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &ab, const BViewType &b);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Pbtrs_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_PBTRS_HPP_
