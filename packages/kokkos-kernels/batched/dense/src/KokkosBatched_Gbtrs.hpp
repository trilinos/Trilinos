// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GBTRS_HPP_
#define KOKKOSBATCHED_GBTRS_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Gbtrs:
///
/// Solve a system of linear equations
///   A * X = B or A**T * X = B or A**H * X = B
///   with a general band matrix A using the LU factorization computed by gbtrf.
/// \tparam ArgTrans: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose) or A**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Gbtrs::Blocked) or unblocked
/// (KokkosBatched::Algo::Gbtrs::Unblocked) algorithm to be used
///
/// \tparam AViewType: Input type for the matrix, needs to be a 2D view
/// \tparam PivViewType: Integer type for pivot indices, needs to be a 1D view
/// \tparam BViewType: Input type for the right-hand side and the solution, needs to be a 1D view
///
/// \param A [in]: A is a ldab by n banded matrix.
/// Details of the LU factorization of the band matrix A, as computed by gbtrf. U is stored as an upper triangular band
/// matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and the multipliers used during the factorization are stored
/// in rows KL+KU+2 to 2*KL+KU+1.
/// \param piv [in]: The pivot indices; for 1 <= i <= N, row i of the matrix was interchanged with row piv(i).
/// \param b [inout]: right-hand side and the solution
/// \param kl [in]: kl specifies the number of subdiagonals within the band of A. kl >= 0
/// \param ku [in]: ku specifies the number of superdiagonals within the band of A. ku >= 0
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgTrans, typename ArgAlgo>
struct SerialGbtrs {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::SerialGbtrs: ArgTrans must be a KokkosBlas::Trans.");
  static_assert(std::is_same_v<ArgAlgo, Algo::Gbtrs::Unblocked>, "KokkosBatched::gbtrs: Use Algo::Gbtrs::Unblocked");
  template <typename AViewType, typename PivViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &piv, const BViewType &b, const int kl,
                                           const int ku);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Gbtrs_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_GBTRS_HPP_
