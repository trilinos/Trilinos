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
#ifndef __KOKKOSBATCHED_FRANCIS_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_FRANCIS_SERIAL_INTERNAL_HPP__

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
struct SerialFrancisInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int mbeg, const int mend, const int morg,
                                           /* */ ValueType *HH, const int hs0, const int hs1,
                                           const Kokkos::complex<ValueType> lambda1,
                                           const Kokkos::complex<ValueType> lambda2, const bool is_complex,
                                           /* */ Kokkos::pair<ValueType, ValueType> *GG, const bool request_schur) {
    typedef ValueType value_type;

    const int hs = hs0 + hs1;
    const value_type one(1), zero(0);
    const Kokkos::pair<value_type, value_type> identity(one, zero);

    /// redefine variables
    const int m = mend - mbeg, mrst = morg - mend, mbeg_mult_hs0 = mbeg * hs0;
    value_type *H = HH + hs * mbeg;

    /// initialize Gs
    Kokkos::pair<value_type, value_type> *Gs = NULL;
    if (request_schur) {
      Gs = (Kokkos::pair<value_type, value_type> *)(GG + mbeg * 2);
      for (int i = 0; i < morg; ++i) {
        GG[2 * i]     = identity;
        GG[2 * i + 1] = identity;
      }
    }

    /// Given a strict Hessenberg matrix H (m x m),
    /// it computes a single implicit QR step with a given shift
    /// - it assumes H has zeros on subdiagonal entries (
    /// givens rotation is defined as G = [gamma -sigma
    ///                                    sigma  gamma]
    ///   G' [chi1 chi2]^t = [alpha 0]^T
    /// where G is stored as a pair of gamma and sigma
    Kokkos::pair<value_type, value_type> G[2];

    /// 0. compute 1st double shift vector
    value_type v[3];
    {
      // this needs m>=3
      // v = M e_1 = (H*H - 2 Re(lambda) H + |lambda|^2 I)e_1
      value_type s, t;
      const value_type h00 = H[0 * hs0 + 0 * hs1], h01 = H[0 * hs0 + 1 * hs1], h10 = H[1 * hs0 + 0 * hs1],
                       h11       = H[1 * hs0 + 1 * hs1],
                       /* */ h21 = H[2 * hs0 + 1 * hs1];
      if (is_complex) {
        s = 2 * lambda1.real();
        t = lambda1.real() * lambda1.real() + lambda1.imag() * lambda1.imag();
      } else {
        const value_type val    = H[(m - 1) * hs];
        const auto dist_lambda1 = Kokkos::ArithTraits<value_type>::abs(lambda1.real() - val);
        const auto dist_lambda2 = Kokkos::ArithTraits<value_type>::abs(lambda2.real() - val);
        const value_type lambda = dist_lambda1 < dist_lambda2 ? lambda1.real() : lambda2.real();
        s                       = 2 * lambda;
        t                       = lambda * lambda;
      }
      v[0] = h00 * h00 + h01 * h10 /* H^2 e_1 */ - s * h00 /* 2 Re(lambda) */ + t;
      v[1] = h10 * h00 + h11 * h10 /*         */ - s * h10;
      v[2] = h21 * h10;
    }

    /// 1. compute the first two givens rotations that introduces a bulge
    {
      SerialGivensInternal::invoke(v[0], v[1], &G[0], &v[0]);
      SerialGivensInternal::invoke(v[0], v[2], &G[1], &v[0]);
      // record
      if (request_schur) {
        Gs[0] = G[0];
        Gs[1] = G[1];
      }

      // apply G' from left and right
      G[0].second = -G[0].second;
      G[1].second = -G[1].second;

      const int mm = m < 4 ? m : 4, nn = m;
      value_type *Hs = H - mbeg_mult_hs0;
      SerialApplyLeftRightGivensInternal ::invoke(G[0], G[1], mm + mbeg, nn + mrst, H, H + hs0, H + 2 * hs0, Hs,
                                                  Hs + hs1, Hs + 2 * hs1, hs0, hs1);
    }

    /// 1. chase the bulge

    // partitions used for loop iteration
    Partition2x2<value_type> H_part2x2(hs0, hs1);
    Partition3x3<value_type> H_part3x3(hs0, hs1);

    // initial partition of A where ATL has a zero dimension
    int m_htl = 1;
    H_part2x2.partWithATL(H, m, m, m_htl, m_htl);
    for (; m_htl < (m - 2); ++m_htl) {
      // part 2x2 into 3x3
      H_part3x3.partWithABR(H_part2x2, 1, 1);
      /// -----------------------------------------------------
      value_type *chi1 = H_part3x3.A11 - hs1;
      value_type *chi2 = chi1 + hs0;
      value_type *chi3 = chi2 + hs0;

      SerialGivensInternal::invoke(*chi1, *chi2, &G[0], chi1);
      *chi2 = zero;
      SerialGivensInternal::invoke(*chi1, *chi3, &G[1], chi1);
      *chi3 = zero;
      // record
      if (request_schur) {
        Gs[2 * m_htl]     = G[0];
        Gs[2 * m_htl + 1] = G[1];
      }

      G[0].second = -G[0].second;
      G[1].second = -G[1].second;

      const int nn = m - m_htl, mtmp = m_htl + 4, mm = mtmp < m ? mtmp : m;
      value_type *a1t = H_part3x3.A11;
      value_type *a2t = a1t + hs0;
      value_type *a3t = a2t + hs0;
      value_type *a1  = H_part3x3.A01 - mbeg_mult_hs0;
      value_type *a2  = a1 + hs1;
      value_type *a3  = a2 + hs1;

      SerialApplyLeftRightGivensInternal ::invoke(G[0], G[1], mm + mbeg, nn + mrst, a1t, a2t, a3t, a1, a2, a3, hs0,
                                                  hs1);
      /// -----------------------------------------------------
      H_part2x2.mergeToATL(H_part3x3);
    }

    // last 3x3 block
    {
      // part 2x2 into 3x3
      H_part3x3.partWithABR(H_part2x2, 1, 1);
      /// -----------------------------------------------------
      value_type *chi1 = H_part3x3.A11 - hs1;
      value_type *chi2 = chi1 + hs0;
      SerialGivensInternal::invoke(*chi1, *chi2, &G[0], chi1);
      *chi2             = zero;
      Gs[2 * m_htl]     = G[0];
      Gs[2 * m_htl + 1] = identity;

      G[0].second = -G[0].second;

      const int mm = m, nn = 2;
      value_type *a1t = H_part3x3.A11;
      value_type *a2t = a1t + hs0;
      value_type *a1  = H_part3x3.A01 - mbeg_mult_hs0;
      value_type *a2  = a1 + hs1;
      SerialApplyLeftRightGivensInternal ::invoke(G[0], mm + mbeg, nn + mrst, a1t, a2t, a1, a2, hs0, hs1);

      /// -----------------------------------------------------
      H_part2x2.mergeToATL(H_part3x3);
    }

    return 0;
  }

  template <typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const int mbeg, const int mend, const int morg,
                                                /* */ ValueType *HH, const int hs0, const int hs1,
                                                const Kokkos::complex<ValueType> lambda1,
                                                const Kokkos::complex<ValueType> lambda2, const bool is_complex) {
    return invoke(mbeg, mend, morg, HH, hs0, hs1, lambda1, lambda2, is_complex,
                  (Kokkos::pair<ValueType, ValueType> *)NULL, false);
  }
};

}  // end namespace KokkosBatched

#endif
