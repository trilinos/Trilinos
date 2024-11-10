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
#ifndef __KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_ApplyHouseholder_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///

struct TeamVectorApplyQ_LeftForwardInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const int k,
                                           /* */ ValueType *A, const int as0, const int as1,
                                           /* */ ValueType *t, const int ts,
                                           /* */ ValueType *B, const int bs0, const int bs1,
                                           /* */ ValueType *w) {
    typedef ValueType value_type;

    /// Given a matrix A that includes a series of householder vectors,
    /// it applies a unitary matrix Q to B from left without transpose
    ///   B = Q B = (H0 H1 H2 H3 ... H(k-1)) B
    /// where
    ///   A is m x k (holding H0, H1 ... H(k-1)
    ///   t is k x 1
    ///   B is m x n

    // partitions used for loop iteration
    Partition2x2<value_type> A_part2x2(as0, as1);
    Partition3x3<value_type> A_part3x3(as0, as1);

    Partition2x1<value_type> t_part2x1(ts);
    Partition3x1<value_type> t_part3x1(ts);

    Partition2x1<value_type> B_part2x1(bs0);
    Partition3x1<value_type> B_part3x1(bs0);

    // initial partition of A where ATL has a zero dimension
    A_part2x2.partWithABR(A, m, k, m - k, 0);
    t_part2x1.partWithAB(t, k, 0);
    B_part2x1.partWithAB(B, m, m - k);

    for (int m_A0 = (k - 1); m_A0 >= 0; --m_A0) {
      // part 2x2 into 3x3
      A_part3x3.partWithATL(A_part2x2, 1, 1);
      t_part3x1.partWithAT(t_part2x1, 1);
      value_type *tau = t_part3x1.A1;

      B_part3x1.partWithAT(B_part2x1, 1);
      const int m_A2 = m - m_A0 - 1;
      /// -----------------------------------------------------
      // left apply householder to partitioned B1 and B2
      TeamVectorApplyLeftHouseholderInternal::invoke(member, m_A2, n, tau, A_part3x3.A21, as0, B_part3x1.A1, bs1,
                                                     B_part3x1.A2, bs0, bs1, w);
      member.team_barrier();
      /// -----------------------------------------------------
      A_part2x2.mergeToABR(A_part3x3);
      t_part2x1.mergeToAB(t_part3x1);
      B_part2x1.mergeToAB(B_part3x1);
    }
    return 0;
  }
};

struct TeamVectorApplyQ_LeftBackwardInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const int k,
                                           /* */ ValueType *A, const int as0, const int as1,
                                           /* */ ValueType *t, const int ts,
                                           /* */ ValueType *B, const int bs0, const int bs1,
                                           /* */ ValueType *w) {
    typedef ValueType value_type;

    /// Given a matrix A that includes a series of householder vectors,
    /// it applies a unitary matrix Q to B from left without transpose
    ///   B = Q B = (H0 H1 H2 H3 ... H(k-1)) B
    /// where
    ///   A is m x k (holding H0, H1 ... H(k-1)
    ///   t is k x 1
    ///   B is m x n

    // partitions used for loop iteration
    Partition2x2<value_type> A_part2x2(as0, as1);
    Partition3x3<value_type> A_part3x3(as0, as1);

    Partition2x1<value_type> t_part2x1(ts);
    Partition3x1<value_type> t_part3x1(ts);

    Partition2x1<value_type> B_part2x1(bs0);
    Partition3x1<value_type> B_part3x1(bs0);

    // initial partition of A where ATL has a zero dimension
    A_part2x2.partWithATL(A, m, k, 0, 0);
    t_part2x1.partWithAT(t, k, 0);
    B_part2x1.partWithAT(B, m, 0);

    for (int m_A0 = 0; m_A0 < k; ++m_A0) {
      // part 2x2 into 3x3
      A_part3x3.partWithABR(A_part2x2, 1, 1);
      t_part3x1.partWithAB(t_part2x1, 1);
      value_type *tau = t_part3x1.A1;

      B_part3x1.partWithAB(B_part2x1, 1);
      const int m_A2 = m - m_A0 - 1;
      /// -----------------------------------------------------
      // left apply householder to partitioned B1 and B2
      TeamVectorApplyLeftHouseholderInternal::invoke(member, m_A2, n, tau, A_part3x3.A21, as0, B_part3x1.A1, bs1,
                                                     B_part3x1.A2, bs0, bs1, w);
      member.team_barrier();
      /// -----------------------------------------------------
      A_part2x2.mergeToATL(A_part3x3);
      t_part2x1.mergeToAT(t_part3x1);
      B_part2x1.mergeToAT(B_part3x1);
    }
    return 0;
  }
};

struct TeamVectorApplyQ_RightForwardInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n, const int k,
                                           /* */ ValueType *A, const int as0, const int as1,
                                           /* */ ValueType *t, const int ts,
                                           /* */ ValueType *B, const int bs0, const int bs1,
                                           /* */ ValueType *w) {
    typedef ValueType value_type;

    /// Given a matrix A that includes a series of householder vectors,
    /// it applies a unitary matrix Q to B from left without transpose
    ///   B = B Q = B (H0 H1 H2 H3 ... H(k-1))
    /// where
    ///   A is n x k (holding H0, H1 ... H(k-1)
    ///   t is k x 1
    ///   B is m x n

    // partitions used for loop iteration
    Partition2x2<value_type> A_part2x2(as0, as1);
    Partition3x3<value_type> A_part3x3(as0, as1);

    Partition2x1<value_type> t_part2x1(ts);
    Partition3x1<value_type> t_part3x1(ts);

    Partition1x2<value_type> B_part1x2(bs1);
    Partition1x3<value_type> B_part1x3(bs1);

    // initial partition of A where ATL has a zero dimension
    A_part2x2.partWithATL(A, n, k, 0, 0);
    t_part2x1.partWithAT(t, k, 0);
    B_part1x2.partWithAL(B, n, 0);

    for (int n_A0 = 0; n_A0 < k; ++n_A0) {
      // part 2x2 into 3x3
      A_part3x3.partWithABR(A_part2x2, 1, 1);
      t_part3x1.partWithAB(t_part2x1, 1);
      value_type *tau = t_part3x1.A1;

      B_part1x3.partWithAR(B_part1x2, 1);
      const int n_B2 = n - n_A0 - 1;
      /// -----------------------------------------------------
      // right apply householder to partitioned B1 and B2
      TeamVectorApplyRightHouseholderInternal::invoke(member, m, n_B2, tau, A_part3x3.A21, as0, B_part1x3.A1, bs0,
                                                      B_part1x3.A2, bs0, bs1, w);
      member.team_barrier();
      /// -----------------------------------------------------
      A_part2x2.mergeToATL(A_part3x3);
      t_part2x1.mergeToAT(t_part3x1);
      B_part1x2.mergeToAL(B_part1x3);
    }
    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
