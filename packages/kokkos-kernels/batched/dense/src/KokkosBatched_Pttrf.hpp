// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PTTRF_HPP_
#define KOKKOSBATCHED_PTTRF_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Pttrf:
/// Compute the Cholesky factorization L*D*L**T (or L*D*L**H) of a real
/// symmetric (or complex Hermitian) positive definite tridiagonal matrix A_l
/// for all l = 0, ..., N
///
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Pttrf::Blocked) or unblocked
/// (KokkosBatched::Algo::Pttrf::Unblocked) algorithm to be used
///
/// \tparam DViewType: Input type for the a diagonal matrix, needs to be a 1D
/// view
/// \tparam EViewType: Input type for the a upper/lower diagonal matrix,
/// needs to be a 1D view
///
/// \param d [inout]: n diagonal elements of the diagonal matrix D
/// \param e [inout]: n-1 upper/lower diagonal elements of the diagonal matrix E
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgAlgo>
struct SerialPttrf {
  static_assert(std::is_same_v<ArgAlgo, Algo::Pttrf::Unblocked>, "KokkosBatched::pttrf: Use Algo::Pttrf::Unblocked");
  template <typename DViewType, typename EViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DViewType &d, const EViewType &e);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Pttrf_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_PTTRF_HPP_
