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
#ifndef __KOKKOSBATCHED_HESSENBERG_QR_WITH_SHIFT_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_HESSENBERG_QR_WITH_SHIFT_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Givens_Serial_Internal.hpp"
#include "KokkosBatched_ApplyGivens_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialHessenbergQR_WithShiftInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int mbeg, const int mend, const int morg,
                                           /* */ ValueType *HH, const int hs0, const int hs1, const ValueType shift,
                                           /* */ Kokkos::pair<ValueType, ValueType> *GG, const bool request_schur) {
    typedef ValueType value_type;
    // typedef Kokkos::ArithTraits<value_type> ats;

    const int hs = hs0 + hs1;
    const value_type zero(0), one(1);
    const Kokkos::pair<value_type, value_type> identity(one, zero);

    /// redefine variables
    const int m = mend - mbeg, mbeg_mult_hs0 = mbeg * hs0;
    value_type *H = HH + mbeg * hs;

    /// initialize Gs
    Kokkos::pair<value_type, value_type> *Gs = NULL;
    if (request_schur) {
      for (int i = 0; i < morg; ++i) GG[i] = identity;
      Gs = GG + mbeg;
    }

    /// Given a strict Hessenberg matrix H (m x m),
    /// it computes a single implicit QR step with a given shift
    /// - it assumes H has zeros on subdiagonal entries (
    /// givens rotation is defined as G = [gamma -sigma
    ///                                    sigma  gamma]
    ///   G' [chi1 chi2]^t = [alpha 0]^T
    /// where G is stored as a pair of gamma and sigma
    Kokkos::pair<value_type, value_type> G;

    /// 0. compute the first givens rotation that introduces a bulge
    {
      const value_type chi1 = H[0] - shift;
      const value_type chi2 = H[hs0];
      /* */ value_type alpha;
      SerialGivensInternal::invoke(chi1, chi2, &G, &alpha);
      // record G
      if (request_schur) Gs[0] = G;

      value_type *h11 = H;
      value_type *h21 = H + hs0;
      value_type *h12 = H + hs1;

      // apply G' from left
      G.second     = -G.second;  // transpose G
      const int nn = m;
      SerialApplyLeftGivensInternal::invoke(G, nn + (morg - mend), h11, hs1, h21, hs1);

      // apply (G')' from right
      const int mm = m < 3 ? m : 3;
      SerialApplyRightGivensInternal::invoke(G, mm + mbeg, h11 - mbeg_mult_hs0, hs0, h12 - mbeg_mult_hs0, hs0);
    }

    /// 1. chase the bulge

    // partitions used for loop iteration
    Partition2x2<value_type> H_part2x2(hs0, hs1);
    Partition3x3<value_type> H_part3x3(hs0, hs1);

    // initial partition of A where ATL has a zero dimension
    H_part2x2.partWithATL(H, m, m, 1, 1);

    for (int m_htl = 1; m_htl < (m - 1); ++m_htl) {
      // part 2x2 into 3x3
      H_part3x3.partWithABR(H_part2x2, 1, 1);
      // const int n_hbr = m - m_htl;
      /// -----------------------------------------------------
      value_type *chi1 = H_part3x3.A11 - hs1;
      value_type *chi2 = H_part3x3.A21 - hs1;
      SerialGivensInternal::invoke(*chi1, *chi2, &G, chi1);
      *chi2 = zero;
      // record G
      if (request_schur) Gs[m_htl] = G;

      G.second = -G.second;  // transpose G

      const int nn = m - m_htl;
      SerialApplyLeftGivensInternal::invoke(G, nn + (morg - mend), H_part3x3.A11, hs1, H_part3x3.A21, hs1);

      const int mtmp = m_htl + 3, mm = mtmp < m ? mtmp : m;
      SerialApplyRightGivensInternal::invoke(G, mm + mbeg, H_part3x3.A01 - mbeg_mult_hs0, hs0,
                                             H_part3x3.A02 - mbeg_mult_hs0, hs0);
      /// -----------------------------------------------------
      H_part2x2.mergeToATL(H_part3x3);
    }
    return 0;
  }

  template <typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const int mbeg, const int mend, const int morg,
                                                /* */ ValueType *HH, const int hs0, const int hs1,
                                                const ValueType shift) {
    return invoke(mbeg, mend, morg, HH, hs0, hs1, shift, (Kokkos::pair<ValueType, ValueType> *)NULL, false);
  }
};

}  // end namespace KokkosBatched

#endif
