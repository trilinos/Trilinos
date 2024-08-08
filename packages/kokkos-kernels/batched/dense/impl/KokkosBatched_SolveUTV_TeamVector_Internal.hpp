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
#ifndef __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas2_team_gemv_impl.hpp"
#include "KokkosBatched_Trsv_TeamVector_Internal.hpp"

#include "KokkosBatched_Gemm_TeamVector_Internal.hpp"
#include "KokkosBatched_Trsm_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Internal
/// ===================
struct TeamVectorSolveUTV_Internal {
  template <typename MemberType, typename ValueType, typename IntType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int matrix_rank, const int m,
                                           const int /*n*/, const ValueType *U, const int us0, const int us1,
                                           const ValueType *T, const int ts0, const int ts1, const ValueType *V,
                                           const int vs0, const int vs1, const IntType *p, const int ps0,
                                           /* */ ValueType *x, const int xs0,
                                           /* */ ValueType *b, const int bs0,
                                           /* */ ValueType *w) {
    typedef ValueType value_type;
    // typedef IntType int_type;

    const value_type one(1), zero(0);
    const int ws0 = 1;

    if (matrix_rank < m) {
      /// w = U^T b
      KokkosBlas::Impl::TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(member, matrix_rank, m, one, U, us1, us0,
                                                                              b, bs0, zero, w, ws0);

      /// w = T^{-1} w
      TeamVectorTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(member, false, matrix_rank, one, T, ts0, ts1, w, ws0);

      /// x = V^T w
      KokkosBlas::Impl::TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(member, m, matrix_rank, one, V, vs1, vs0,
                                                                              w, ws0, zero, x, xs0);
    } else {
      KokkosBlas::Impl::TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(member, matrix_rank, m, one, U, us1, us0,
                                                                              b, bs0, zero, x, xs0);

      TeamVectorTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(member, false, matrix_rank, one, T, ts0, ts1, x, xs0);
    }

    /// x = P^T x
    TeamVectorApplyPivotVectorBackwardInternal ::invoke(member, m, p, ps0, x, xs0);

    return 0;
  }

  template <typename MemberType, typename ValueType, typename IntType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int matrix_rank, const int m, const int n,
                                           const int nrhs, const ValueType *U, const int us0, const int us1,
                                           const ValueType *T, const int ts0, const int ts1, const ValueType *V,
                                           const int vs0, const int vs1, const IntType *p, const int ps0,
                                           /* */ ValueType *X, const int xs0, const int xs1,
                                           /* */ ValueType *B, const int bs0, const int bs1,
                                           /* */ ValueType *w) {
    typedef ValueType value_type;
    // typedef IntType int_type;

    const value_type one(1), zero(0);

    value_type *W = w;  /// m x nrhs
    const int ws0 = xs0 < xs1 ? 1 : nrhs, ws1 = xs0 < xs1 ? m : 1;

    if (matrix_rank < n) {
      /// U is m x matrix_rank
      /// T is matrix_rank x matrix_rank
      /// V is matrix_rank x n
      /// W = U^T B
      TeamVectorGemmInternal<Algo::Gemm::Unblocked>::invoke(member, matrix_rank, nrhs, m, one, U, us1, us0, B, bs0, bs1,
                                                            zero, W, ws0, ws1);
      member.team_barrier();

      /// W = T^{-1} W
      TeamVectorTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(member, false, matrix_rank, nrhs, one, T, ts0, ts1,
                                                                     W, ws0, ws1);
      member.team_barrier();

      /// X = V^T W
      TeamVectorGemmInternal<Algo::Gemm::Unblocked>::invoke(member, n, nrhs, matrix_rank, one, V, vs1, vs0, W, ws0, ws1,
                                                            zero, X, xs0, xs1);
      member.team_barrier();
    } else {
      /// W = U^T B
      TeamVectorGemmInternal<Algo::Gemm::Unblocked>::invoke(member, matrix_rank, nrhs, m, one, U, us1, us0, B, bs0, bs1,
                                                            zero, X, xs0, xs1);
      member.team_barrier();

      /// X = T^{-1} X
      TeamVectorTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, false, matrix_rank, nrhs, one, T, ts0, ts1,
                                                                     X, xs0, xs1);
      member.team_barrier();
    }

    /// X = P^T X
    TeamVectorApplyPivotMatrixBackwardInternal ::invoke(member, nrhs, n, p, ps0, X, xs0, xs1);

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
