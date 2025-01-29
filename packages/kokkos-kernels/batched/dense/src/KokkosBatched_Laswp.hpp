//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef KOKKOSBATCHED_LASWP_HPP_
#define KOKKOSBATCHED_LASWP_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Laswp:
///
/// performs a series of row interchanges on the matrix A.
/// One row interchange is initiated for each of rows K1 through K2 of A.
///
/// \tparam PivViewType: Input type for the a superdiagonal matrix, needs to
/// be a 1D view
/// \tparam AViewType: Input type for the vector or matrix, needs to be a 1D or
/// 2D view
///
/// \param piv [in]: The pivot indices; for 0 <= i < N, row i of the
/// matrix was interchanged with row piv(i).
/// \param A [inout]: A is a lda by n matrix. The matrix of column dimension N
/// to which the row interchanges will be applied.
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgDirect>
struct SerialLaswp {
  template <typename PivViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const PivViewType &piv, const AViewType &A);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Laswp_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_LASWP_HPP_
