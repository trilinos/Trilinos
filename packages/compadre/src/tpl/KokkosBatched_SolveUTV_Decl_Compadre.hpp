// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef __KOKKOSBATCHED_SOLVE_UTV_DECL_COMPADRE_HPP__
#define __KOKKOSBATCHED_SOLVE_UTV_DECL_COMPADRE_HPP__

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
  /// When A is a full rank i.e., matrix_rank == m, UTV computes QR with column pivoting only
  /// where Q is stored in U and R is stored in T
  ///
  
  ///
  /// TeamVector Solve UTV
  ///

  template<typename MemberType,
           typename ArgAlgo>
  struct TeamVectorSolveUTVCompadre {
    template<typename UViewType,
         typename TViewType,
         typename VViewType,
         typename pViewType,
         typename BViewType,
         typename XViewType,
         typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
       const int matrix_rank,
           const UViewType &U,
       const TViewType &T,
       const VViewType &V,
       const pViewType &p,
       const BViewType &B,
       const XViewType &X,
       const wViewType &w_a,
       const wViewType &w_b,
       const bool implicit_RHS);
  };

}

#include "KokkosBatched_SolveUTV_TeamVector_Impl_Compadre.hpp"

#endif
