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
/// \tparam ABViewType: Input type for a banded matrix, needs to be a 2D
/// view
/// \tparam BViewType: Input type for a right-hand side and the solution,
/// needs to be a 1D view
///
/// \param ab [in]: ab is a ldab by n banded matrix, with ( kd + 1 ) diagonals
/// \param b  [inout]: right-hand side and the solution, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgUplo, typename ArgAlgo>
struct SerialPbtrs {
  template <typename ABViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &ab, const BViewType &b);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Pbtrs_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_PBTRS_HPP_
