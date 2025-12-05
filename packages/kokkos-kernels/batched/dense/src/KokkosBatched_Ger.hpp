// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GER_HPP_
#define KOKKOSBATCHED_GER_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Ger:
/// Performs the rank 1 operation
///   A := alpha*x*y**T + A or A := alpha*x*y**H + A
///    where alpha is a scalar, x is an m element vector, y is an n element
/// vector and A is an m by n matrix.
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \tparam YViewType: Input type for the vector y, needs to be a 1D view
/// \tparam AViewType: Input/output type for the matrix A, needs to be a 2D view
///
/// \param alpha [in]: A is a m by n general matrix, a rank 2 view
/// \param x [in]: x is a length m vector, a rank 1 view
/// \param y [in]: y is a length n vector, a rank 1 view
/// \param A [inout]: A is a m by n matrix, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgTrans>
struct SerialGer {
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &a);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Ger_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_GER_HPP_
