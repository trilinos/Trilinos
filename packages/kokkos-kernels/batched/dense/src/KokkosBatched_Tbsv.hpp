// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_TBSV_HPP_
#define KOKKOSBATCHED_TBSV_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Tbsv:
///
/// Solve Ab_l x_l = b_l for all l = 0, ..., N
///   using the triangular solve algorithm Tbsv. Ab is an n by n unit, or
///   non-unit, upper or lower triangular band matrix, with ( k + 1 )
///   diagonals.
///
/// \tparam ArgUplo: Type indicating whether A is the upper (Uplo::Upper) or lower (Uplo::Lower) triangular matrix
/// \tparam ArgTrans: Type indicating the equations to be solved as follows
///  - ArgTrans::NoTranspose:   A * X = B
///  - ArgTrans::Transpose:     A**T * X = B
///  - ArgTrans::ConjTranspose: A**H * X = B
/// \tparam ArgDiag: Type indicating whether A is the unit (Diag::Unit) or non-unit (Diag::NonUnit) triangular matrix
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Tbsv::Blocked) or unblocked
/// (KokkosBatched::Algo::Tbsv::Unblocked) algorithm to be used
///
/// \tparam AViewType: Input type for the matrix, needs to be a 2D view
/// \tparam XViewType: Input type for the right-hand side and the solution,
/// needs to be a 1D view
///
/// \param A [in]: A is a lda by n banded matrix, with ( k + 1 ) diagonals
/// \param X [inout]: right-hand side and the solution, a rank 1 view
/// \param k [in]: k specifies the number of superdiagonals or subdiagonals of
/// matrix A. k >= 0
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
struct SerialTbsv {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::tbsv: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::is_same_v<ArgTrans, Trans::NoTranspose> || std::is_same_v<ArgTrans, Trans::Transpose> ||
                    std::is_same_v<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::tbsv: Use Trans::NoTranspose, Trans::Transpose or Trans::ConjTranspose");
  static_assert(
      std::is_same_v<ArgDiag, Diag::Unit> || std::is_same_v<ArgDiag, Diag::NonUnit>,
      "KokkosBatched::tbsv: Use Diag::Unit for unit triangular matrix or Diag::NonUnit for non-unit triangular matrix");
  static_assert(std::is_same_v<ArgAlgo, Algo::Tbsv::Unblocked>, "KokkosBatched::tbsv: Use Algo::Tbsv::Unblocked");
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &X, const int k);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Tbsv_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_TBSV_HPP_
