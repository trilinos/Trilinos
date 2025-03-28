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
#ifndef KOKKOSBATCHED_GETRS_HPP_
#define KOKKOSBATCHED_GETRS_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Getrs:
/// Solve a system of linear equations
///   A * x = b or A**T * x = b
///   with a general N-by-N matrix A using LU factorization computed
///   by Getrf.
/// \tparam AViewType: Input type for the matrix, needs to be a 2D view
/// \tparam PivViewType: Input type for the pivot indices, needs to be a 1D view
/// \tparam BViewType: Input type for the right-hand side and the solution,
/// needs to be a 1D view
///
/// \param A [in]: A is a m by n general matrix, a rank 2 view
/// \param piv [in]: On exit, the pivot indices, a rank 1 view
/// \param B [inout]: right-hand side and the solution, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgTrans, typename ArgAlgo>
struct SerialGetrs {
  template <typename AViewType, typename PivViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &piv, const BViewType &b);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Getrs_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_GETRF_HPP_
