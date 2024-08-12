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
#ifndef __KOKKOSBATCHED_LEFT_EIGENVECTOR_FROM_SCHUR_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_LEFT_EIGENVECTOR_FROM_SCHUR_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_ShiftedTrsv_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialLeftEigenvectorFromSchurInternal {
  /// Given a quasi upper triangular matrix S (m x m), this computes all left
  /// eigenvectors.
  ///
  /// Parameters:
  ///   [in]m
  ///     A dimension of the square matrix S.
  ///   [in]S, [in]ss0, [in]ss1
  ///     A quasi upper triangular part of Schur decomposition which is computed
  ///       A = U^H S U
  ///   [out]V, [in]vs0, [out]vs1
  ///     A set of left eigen vectors.
  ///   [in]w
  ///     contiguous workspace that can hold complex array (m)
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m,
                                           /* */ ValueType *S, const int ss0, const int ss1,
                                           /* */ ValueType *V, const int vs0, const int vs1,
                                           /* */ ValueType *w, const int *blks) {
    typedef ValueType value_type;
    typedef Kokkos::ArithTraits<value_type> ats;
    // typedef typename ats::mag_type mag_type;
    typedef Kokkos::complex<value_type> complex_type;

    const value_type zero(0), one(1);
    /// SerialSetInternal::invoke(m, m, zero, V, vs0, vs1);

    value_type *b = w;  // consider complex case

    /// partitions used for loop iteration
    Partition2x2<value_type> S_part2x2(ss0, ss1);
    Partition3x3<value_type> S_part3x3(ss0, ss1);

    Partition2x1<value_type> V_part2x1(vs0);
    Partition3x1<value_type> V_part3x1(vs0);

    /// initial partition of S where ATL has a zero dimension
    S_part2x2.partWithATL(S, m, m, 0, 0);
    V_part2x1.partWithAT(V, m, 0);

    // const mag_type tol = ats::epsilon();
    int m_stl = 0;
    for (; m_stl < (m - 1);) {
      /// part 2x2 into 3x3
      const int mA11 = blks[m_stl];
      assert(((mA11 == 1) || (mA11 == 2)) && "LeftEigenvectorFromSchur: blk is not 1x1 nor 2x2");

      S_part3x3.partWithABR(S_part2x2, mA11, mA11);
      V_part3x1.partWithAB(V_part2x1, mA11);

      const int m_stl_plus_mA11 = m_stl + mA11;
      if (mA11 == 1) {
        /// real eigenvalue
        const value_type lambda = *S_part3x3.A11;

        /// initialize a left hand side
        b[m_stl] = one;
        for (int j = 0; j < (m - m_stl_plus_mA11); ++j) b[j + m_stl_plus_mA11] = -S_part3x3.A12[j * ss1];

        /// perform shifted trsv (transposed)
        SerialShiftedTrsvInternalLower::invoke(m - m_stl_plus_mA11, lambda, S_part3x3.A22, ss1, ss0,
                                               b + m_stl_plus_mA11, 1, blks + m_stl_plus_mA11);

        /// copy back to V (row wise copy)
        for (int j = 0; j < m_stl; ++j) V_part3x1.A1[j * vs1] = zero;
        for (int j = m_stl; j < m; ++j) V_part3x1.A1[j * vs1] = b[j];
      } else {
        /// complex eigen pair
        const value_type alpha11 = S_part3x3.A11[0], alpha12 = S_part3x3.A11[ss1], alpha21 = S_part3x3.A11[ss0],
                         beta = ats::sqrt(-alpha12 * alpha21);

        const complex_type lambda(alpha11, beta);
        complex_type *bc = (complex_type *)(b);

        /// initialize a left hand side
        bc[m_stl]     = complex_type(beta, zero);
        bc[m_stl + 1] = complex_type(zero, -alpha12);

        const value_type *S_A12_a = S_part3x3.A12;
        const value_type *S_A12_b = S_part3x3.A12 + ss0;
        for (int j = 0; j < (m - m_stl_plus_mA11); ++j)
          bc[j + m_stl_plus_mA11] = complex_type(-S_A12_a[j * ss1] * beta, S_A12_b[j * ss1] * alpha12);

        /// perform shifted trsv
        SerialShiftedTrsvInternalLower::invoke(m - m_stl_plus_mA11, lambda, S_part3x3.A22, ss1, ss0,
                                               bc + m_stl_plus_mA11, 1, blks + m_stl_plus_mA11);

        /// copy back to V
        value_type *V_A1_r = V_part3x1.A1;
        value_type *V_A1_i = V_part3x1.A1 + vs0;
        for (int j = 0; j < m_stl; ++j) {
          V_A1_r[j * vs1] = zero;
          V_A1_i[j * vs1] = zero;
        }
        for (int j = m_stl; j < m; ++j) {
          V_A1_r[j * vs1] = bc[j].real();
          V_A1_i[j * vs1] = bc[j].imag();
        }
        /// ---------------------------------------------------
      }
      S_part2x2.mergeToATL(S_part3x3);
      V_part2x1.mergeToAT(V_part3x1);
      m_stl += mA11;
    }

    /// case: m_stl = m-1
    if (m_stl < m) {
      value_type *VV = V + m_stl * vs0;
      for (int j = 0; j < m_stl; ++j) VV[j * vs1] = zero;
      VV[m_stl * vs1] = one;
    }

    return 0;
  }
};

}  // namespace KokkosBatched

#endif
