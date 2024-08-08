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
#ifndef __KOKKOSBATCHED_HESSENBERG_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_HESSENBERG_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_Serial_Internal.hpp"
#include "KokkosBatched_ApplyHouseholder_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialHessenbergInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m,  // m = NumRows(A)
                                           const int n,  // n = NumCols(A)
                                           /* */ ValueType *A, const int as0, const int as1,
                                           /* */ ValueType *t, const int ts,
                                           /* */ ValueType *w) {
    typedef ValueType value_type;

    /// Given a matrix A, it computes hessenberg decomposition of the matrix
    ///  - t is to store tau and w is for workspace
    ///  - H = Q^H A Q and A = Q H Q^H

    // partitions used for loop iteration
    Partition2x2<value_type> A_part2x2(as0, as1);
    Partition3x3<value_type> A_part3x3(as0, as1);

    Partition2x1<value_type> t_part2x1(ts);
    Partition3x1<value_type> t_part3x1(ts);

    // partitions used in loop body
    Partition2x1<value_type> A21_part2x1(as0);
    Partition2x1<value_type> A22_part2x1(as0);
    Partition1x2<value_type> A2_part1x2(as1);

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
      if (m_A22 > 0) {
        // partition A21 into 2x1
        A21_part2x1.partWithAT(A_part3x3.A21, m_A22, 1);

        // perform householder transformation
        const int m_A22_b = m_A22 - 1;
        SerialLeftHouseholderInternal::invoke(m_A22_b, A21_part2x1.AT, A21_part2x1.AB, as0, tau);

        // partition A22 into 2x1
        A22_part2x1.partWithAT(A_part3x3.A22, m_A22, 1);

        // left apply householder to partitioned A22
        SerialApplyLeftHouseholderInternal::invoke(m_A22_b, n_A22, tau, A21_part2x1.AB, as0, A22_part2x1.AT, as1,
                                                   A22_part2x1.AB, as0, as1, w);

        // partition A*2 column into 1x2
        A2_part1x2.partWithAL(A_part3x3.A02, n_A22, 1);

        // right apply householder to A*2 colums
        const int n_A22_r = n_A22 - 1;
        SerialApplyRightHouseholderInternal::invoke(m, n_A22_r, tau, A21_part2x1.AB, as0, A2_part1x2.AL, as0,
                                                    A2_part1x2.AR, as0, as1, w);
      }
      /// -----------------------------------------------------
      A_part2x2.mergeToATL(A_part3x3);
      t_part2x1.mergeToAT(t_part3x1);
    }
    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
