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
#ifndef KOKKOSBATCHED_PBTRF_HPP_
#define KOKKOSBATCHED_PBTRF_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Pbtrf:
/// Compute the Cholesky factorization U**H * U (or L * L**H) of a real
/// symmetric (or complex Hermitian) positive definite banded matrix A_l
/// for all l = 0, ...,
/// The factorization has the form
///    A = U**T * U ,  if ArgUplo = KokkosBatched::Uplo::Upper, or
///    A = L  * L**T,  if ArgUplo = KokkosBatched::Uplo::Lower,
/// where U is an upper triangular matrix, U**T is the transpose of U, and
/// L is lower triangular.
/// This is the unblocked version of the algorithm, calling Level 2 BLAS.
///
/// \tparam ArgUplo: Type indicating whether A is the upper (Uplo::Upper) or lower (Uplo::Lower) triangular matrix
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Pbtrf::Blocked) or unblocked
/// (KokkosBatched::Algo::Pbtrf::Unblocked) algorithm to be used
///
/// \tparam ABViewType: Input type for a banded matrix, needs to be a 2D view
///
/// \param ab [inout]: ab is a ldab by n banded matrix, with ( kd + 1 ) diagonals
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgUplo, typename ArgAlgo>
struct SerialPbtrf {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::pbtrf: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::is_same_v<ArgAlgo, Algo::Pbtrf::Unblocked>, "KokkosBatched::pbtrf: Use Algo::Pbtrf::Unblocked");
  template <typename ABViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ABViewType &ab);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Pbtrf_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_PBTRF_HPP_
