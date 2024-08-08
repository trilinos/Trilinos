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
#ifndef __KOKKOSBATCHED_RIGHT_EIGENVECTOR_FROM_SCHUR_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_RIGHT_EIGENVECTOR_FROM_SCHUR_SERIAL_INTERNAL_HPP__

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
struct SerialRightEigenvectorFromSchurInternal {
  /// Given a quasi upper triangular matrix S (m x m), this computes all right
  /// eigenvectors.
  ///
  /// Parameters:
  ///   [in]m
  ///     A dimension of the square matrix S.
  ///   [in]S, [in]ss0, [in]ss1
  ///     A quasi upper triangular part of Schur decomposition which is computed
  ///       A = U^H S U
  ///   [out]V, [in]vs0, [out]vs1
  ///     A set of right eigen vectors.
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
    // const int ss(ss0+ss1);
    /// SerialSetInternal::invoke(m, m, zero, V, vs0, vs1);

    value_type *b = w;  // consider complex case

    /// partitions used for loop iteration
    Partition2x2<value_type> S_part2x2(ss0, ss1);
    Partition3x3<value_type> S_part3x3(ss0, ss1);

    Partition1x2<value_type> V_part1x2(vs1);
    Partition1x3<value_type> V_part1x3(vs1);

    /// initial partition of S where ABR has a zero dimension
    S_part2x2.partWithABR(S, m, m, 0, 0);
    V_part1x2.partWithAR(V, m, 0);

    // const mag_type tol = ats::epsilon();
    int m_stl = m;
    for (; m_stl > 0;) {
      /// part 2x2 into 3x3
      const int mA11 = blks[m_stl - 1];
      assert(((mA11 == 1) || (mA11 == 2)) && "RightEigenvectorFromSchur: blk is not 1x1 nor 2x2");

      S_part3x3.partWithATL(S_part2x2, mA11, mA11);
      V_part1x3.partWithAL(V_part1x2, mA11);

      const int m_stl_minus_mA11 = m_stl - mA11;
      if (mA11 == 1) {
        /// real eigenvalue
        const value_type lambda = *S_part3x3.A11;

        /// initialize a right eigen vector
        for (int i = 0; i < m_stl_minus_mA11; ++i) b[i] = -S_part3x3.A01[i * ss0];
        b[m_stl - 1] = one;

        /// perform shifted trsv
        SerialShiftedTrsvInternalUpper::invoke(m_stl_minus_mA11, lambda, S_part3x3.A00, ss0, ss1, w, 1, blks);

        /// copy back to V
        for (int i = 0; i < m_stl; ++i) V_part1x3.A1[i * vs0] = w[i];
        for (int i = m_stl; i < m; ++i) V_part1x3.A1[i * vs0] = zero;
      } else {
        /// complex eigen pair
        const value_type alpha11 = S_part3x3.A11[0], alpha12 = S_part3x3.A11[ss1], alpha21 = S_part3x3.A11[ss0],
                         beta = ats::sqrt(-alpha12 * alpha21);

        const complex_type lambda(alpha11, beta);
        complex_type *bc = (complex_type *)(b);

        /// initialize a right eigen vector
        const value_type *S_A01_a = S_part3x3.A01;
        const value_type *S_A01_b = S_part3x3.A01 + ss1;
        for (int i = 0; i < m_stl_minus_mA11; ++i)
          bc[i] = complex_type(-S_A01_a[i * ss0] * beta, S_A01_b[i * ss0] * alpha21);
        bc[m_stl - 2] = complex_type(beta, zero);
        bc[m_stl - 1] = complex_type(zero, -alpha21);

        /// perform shifted trsv
        SerialShiftedTrsvInternalUpper::invoke(m_stl_minus_mA11, lambda, S_part3x3.A00, ss0, ss1, bc, 1, blks);

        /// copy back to V
        value_type *V_A1_r = V_part1x3.A1;
        value_type *V_A1_i = V_part1x3.A1 + vs1;
        for (int i = 0; i < m_stl; ++i) {
          V_A1_r[i * vs0] = bc[i].real();
          V_A1_i[i * vs0] = bc[i].imag();
        }
        for (int i = m_stl; i < m; ++i) {
          V_A1_r[i * vs0] = zero;
          V_A1_i[i * vs0] = zero;
        }
        /// ---------------------------------------------------
      }
      S_part2x2.mergeToABR(S_part3x3);
      V_part1x2.mergeToAR(V_part1x3);
      m_stl -= mA11;
    }
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
