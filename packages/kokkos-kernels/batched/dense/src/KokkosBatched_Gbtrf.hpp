// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GBTRF_HPP_
#define KOKKOSBATCHED_GBTRF_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Gbtrf:
/// Compute a LU factorization of a m-by-n band matrix A using partial
///   pivoting with row interchanges.
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Gbtrf::Blocked) or unblocked
/// (KokkosBatched::Algo::Gbtrf::Unblocked) algorithm to be used
///
/// \tparam ABViewType: Input type for the matrix, needs to be a 2D view
/// \tparam PivViewType: Input type for the pivot indices, needs to be a 1D view
///
/// \param Ab [inout]: A is a ldab by n band matrix, a rank 2 view
/// \param piv [out]: On exit, the pivot indices, a rank 1 view
/// \param kl [in]: The number of subdiagonals within the band of A. kl >= 0
/// \param ku [in]: The number of superdiagonals within the band of A. ku >= 0
/// \param m [in]: The number of rows of the matrix A. (optional, default is -1, corresponding to m == n)
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgAlgo>
struct SerialGbtrf {
  static_assert(std::is_same_v<ArgAlgo, Algo::Gbtrf::Unblocked>, "KokkosBatched::gbtrf: Use Algo::Gbtrf::Unblocked");
  template <typename ABViewType, typename PivViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &Ab, const PivViewType &piv, const int kl, const int ku,
                                           const int m = -1);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Gbtrf_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_GBTRF_HPP_
