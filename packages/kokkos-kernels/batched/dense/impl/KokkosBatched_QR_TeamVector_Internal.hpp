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
#ifndef __KOKKOSBATCHED_QR_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_QR_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_TeamVector_Internal.hpp"
#include "KokkosBatched_ApplyHouseholder_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct TeamVectorQR_Internal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const int m,  // m = NumRows(A)
                                           const int n,  // n = NumCols(A)
                                           /* */ ValueType *A, const int as0, const int as1,
                                           /* */ ValueType *t, const int ts,
                                           /* */ ValueType *w) {
    typedef ValueType value_type;

    /// Given a matrix A, it computes QR decomposition of the matrix
    ///  - t is to store tau and w is for workspace

    // partitions used for loop iteration
    Partition2x2<value_type> A_part2x2(as0, as1);
    Partition3x3<value_type> A_part3x3(as0, as1);

    Partition2x1<value_type> t_part2x1(ts);
    Partition3x1<value_type> t_part3x1(ts);

    // initial partition of A where ATL has a zero dimension
    A_part2x2.partWithATL(A, m, n, 0, 0);
    t_part2x1.partWithAT(t, m, 0);

    for (int m_atl = 0; m_atl < m; ++m_atl) {
      // part 2x2 into 3x3
      A_part3x3.partWithABR(A_part2x2, 1, 1);
      const int m_A22 = m - m_atl - 1;
      const int n_A22 = n - m_atl - 1;

      t_part3x1.partWithAB(t_part2x1, 1);
      value_type *tau = t_part3x1.A1;

      /// -----------------------------------------------------

      // perform householder transformation
      TeamVectorLeftHouseholderInternal::invoke(member, m_A22, A_part3x3.A11, A_part3x3.A21, as0, tau);
      member.team_barrier();

      // left apply householder to A22
      TeamVectorApplyLeftHouseholderInternal::invoke(member, m_A22, n_A22, tau, A_part3x3.A21, as0, A_part3x3.A12, as1,
                                                     A_part3x3.A22, as0, as1, w);
      member.team_barrier();
      /// -----------------------------------------------------
      A_part2x2.mergeToATL(A_part3x3);
      t_part2x1.mergeToAT(t_part3x1);
    }
    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
