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
#ifndef __KOKKOSBATCHED_UTV_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_UTV_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_SetTriangular_Internal.hpp"
#include "KokkosBatched_QR_TeamVector_Internal.hpp"
#include "KokkosBatched_QR_WithColumnPivoting_TeamVector_Internal.hpp"
#include "KokkosBatched_QR_FormQ_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Internal
/// ===================
struct TeamVectorUTV_Internal {
  template <typename MemberType, typename ValueType, typename IntType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m,
                                           const int n,  // m = NumRows(A), n = NumCols(A)
                                           /* */ ValueType *A, const int as0, const int as1,
                                           /* */ IntType *p, const int ps0,
                                           /* */ ValueType *U, const int us0, const int us1,
                                           /* */ ValueType *V, const int vs0, const int vs1,
                                           /* */ ValueType *w,  // 3*m, tau, norm, householder workspace
                                           /* */ int &matrix_rank) {
    typedef ValueType value_type;
    // typedef IntType int_type;

    value_type *t = w;
    w += m;
    const int ts0(1);

    value_type *work = w;

    matrix_rank = -1;
    TeamVectorQR_WithColumnPivotingInternal ::invoke(member, m, n, A, as0, as1, t, ts0, p, ps0, work, matrix_rank);

    TeamVectorQR_FormQ_Internal ::invoke(member, m, matrix_rank, matrix_rank, A, as0, as1, t, ts0, U, us0, us1, work);
    member.team_barrier();

    /// for rank deficient matrix
    if (matrix_rank < n) {
      const value_type zero(0);
      TeamVectorSetLowerTriangularInternal ::invoke(member, matrix_rank, matrix_rank, 1, zero, A, as0, as1);

      TeamVectorQR_Internal ::invoke(member, n, matrix_rank, A, as1, as0, t, ts0, work);

      TeamVectorQR_FormQ_Internal ::invoke(member, n, matrix_rank, matrix_rank, A, as1, as0, t, ts0, V, vs1, vs0, work);
    }

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
