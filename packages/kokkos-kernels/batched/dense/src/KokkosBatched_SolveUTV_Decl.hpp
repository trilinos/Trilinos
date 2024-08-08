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
#ifndef __KOKKOSBATCHED_SOLVE_UTV_DECL_HPP__
#define __KOKKOSBATCHED_SOLVE_UTV_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// For given UTV = A P^T, it solves A X = B
/// - input:
///   - matrix_rank is computed while UTV factorization
///   - U is m x m real matrix (m x matrix_rank is only used)
///   - T is m x m real matrix (matrix_rank x matrix_rank is only used)
///   - V is m x m real matrix (matrix_Rank x m is only used)
///   - p is m integer vector including pivot indicies
///   - X is a solution matrix (or vector)
///   - B is a right hand side matrix (or vector)
///   - w is B.span() real vector workspace (contiguous)
/// - output:
///   - B is overwritten with its solutions
///
/// When A is a full rank i.e., matrix_rank == m, UTV computes QR with column
/// pivoting only where Q is stored in U and R is stored in T
///

///
/// TeamVector Solve UTV
///

template <typename MemberType, typename ArgAlgo>
struct TeamVectorSolveUTV {
  template <typename UViewType, typename TViewType, typename VViewType, typename pViewType, typename XViewType,
            typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int matrix_rank, const UViewType &U,
                                           const TViewType &T, const VViewType &V, const pViewType &p,
                                           const XViewType &X, const BViewType &B, const wViewType &w);
};

}  // namespace KokkosBatched

#include "KokkosBatched_SolveUTV_TeamVector_Impl.hpp"

#endif
