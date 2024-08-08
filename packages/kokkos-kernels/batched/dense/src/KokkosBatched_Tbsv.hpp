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
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &X, const int k);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Tbsv_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_TBSV_HPP_
