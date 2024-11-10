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
#ifndef __KOKKOSBATCHED_EIGENVALUE_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_EIGENVALUE_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_WilkinsonShift_Serial_Internal.hpp"
#include "KokkosBatched_Schur2x2_Serial_Internal.hpp"
#include "KokkosBatched_HessenbergQR_WithShift_Serial_Internal.hpp"
#include "KokkosBatched_Francis_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialEigenvalueInternal {
  /// Given a strictly Hessenberg matrix H (m x m), this computes all
  /// eigenvalues using the Francis method and stores them into a vector e. This
  /// routine does not scale nor balance the matrix for the numerical stability.
  ///
  /// Parameters:
  ///   [in]m
  ///     A dimension of the square matrix H.
  ///   [in/out]H, [in]hs0, [in]hs1
  ///     Real Hessenberg matrix H(m x m) with strides hs0 and hs1.
  ///     Entering the routine, H is assumed to have a upper Hessenberg form,
  ///     where all subdiagonals are zero. The matrix is overwritten on exit.
  ///   [out]er, [in]ers, [out]ei, [in]eis
  ///     A complex vector er(m)+ei(m)i with a stride ers and eis to store
  ///     computed eigenvalues. For a complex eigen pair, it stores a+bi and
  ///     a-bi consecutively.
  ///   [in]restart(false)
  ///     With a restart option, the routine assume that the matrix H and the
  ///     vector e contain the partial results from the previous run. When m = 1
  ///     or 2, this option won't work as the routine always computes the all
  ///     eigenvalues.
  ///   [in]max_iteration(300)
  ///     Unlike LAPACK which uses various methods for different types of
  ///     matrices, this routine uses the Francis method only. A user can set
  ///     the maximum number of iterations. When it reaches the maximum
  ///     iteration counts without converging all eigenvalues, the routine
  ///     returns -1.
  template <typename RealType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m,
                                           /* */ RealType *H, const int hs0, const int hs1,
                                           /* */ RealType *er, const int ers,
                                           /* */ RealType *ei, const int eis, const bool restart = false,
                                           const int user_max_iteration = -1) {
    typedef RealType real_type;
    typedef Kokkos::ArithTraits<real_type> ats;
    const real_type zero(0), nan(ats::nan()), tol = 1e2 * ats::epsilon();
    const int max_iteration = user_max_iteration < 0 ? 300 : user_max_iteration;

    int r_val = 0;
    if (restart) {
      if (m <= 2) {
        Kokkos::abort("Error: restart option cannot be used for m=1 or m=2");
      }
    } else {
      for (int i = 0; i < m; ++i) er[i * ers] = nan;
    }

    const int hs = hs0 + hs1;  /// diagonal stride
    switch (m) {
      case 0: { /* do nothing */ break;
      }
      case 1: {
        er[0] = H[0];
        ei[0] = zero;
        break;
      }
      case 2: {
        /// compute eigenvalues from the characteristic determinant equation
        bool is_complex;
        Kokkos::complex<real_type> lambda1, lambda2;
        SerialWilkinsonShiftInternal::invoke(H[0], H[hs1], H[hs0], H[hs], &lambda1, &lambda2, &is_complex);
        er[0] = lambda1.real();
        ei[0] = lambda1.imag();
        er[1] = lambda2.real();
        ei[1] = lambda2.imag();
        break;
      }
      default: {
        /// Francis method
        int iter(0);            /// iteration count
        bool converge = false;  /// bool to check all eigenvalues are converged

        while (!converge && iter < max_iteration) {
          /// Step 1: find a set of unrevealed eigenvalues
          int cnt = 1;

          /// find mbeg (first nonzero subdiag value)
          for (; cnt < m; ++cnt) {
            const auto val = ats::abs(*(H + cnt * hs - hs1));
            if (val > tol) break;
          }
          const int mbeg = cnt - 1;

          /// find mend (first zero subdiag value)
          for (; cnt < m; ++cnt) {
            const auto val = ats::abs(*(H + cnt * hs - hs1));
            if (val < tol) break;
          }
          const int mend  = cnt;
          const int mdiff = mend - mbeg;

          /// Step 2: if there exist non-converged eigen values
          if (1 < mdiff) {
#if 0  /// implicit QR with shift for testing 
            {
              /// Rayleigh quotient shift 
              const real_type shift = *(H+(mend-1)*hs); 
              SerialHessenbergQR_WithShiftInternal::invoke(0, mdiff, mdiff, 
                                                           H+hs*mbeg, hs0, hs1,
                                                           shift);
              real_type *sub2x2 = H+(mend-2)*hs;
              const auto     val = ats::abs(sub2x2[hs0]);
              if (val < tol) { /// this eigenvalue converges
                er[(mend-1)*ers] = sub2x2[hs]; ei[(mend-1)*eis] = zero;
              }
            }
#endif

#if 1  /// Francis double shift method
            {
              /// find a complex eigen pair
              Kokkos::complex<real_type> lambda1, lambda2;
              bool is_complex;
              real_type *sub2x2 = H + (mend - 2) * hs;
              if (2 == mdiff) {
                SerialWilkinsonShiftInternal::invoke(sub2x2[0], sub2x2[hs1], sub2x2[hs0], sub2x2[hs], &lambda1,
                                                     &lambda2, &is_complex);
                sub2x2[hs0] = zero;

                /// eigenvalues are from wilkinson shift
                er[(mbeg + 0) * ers] = lambda1.real();
                ei[(mbeg + 0) * eis] = lambda1.imag();
                er[(mbeg + 1) * ers] = lambda2.real();
                ei[(mbeg + 1) * eis] = lambda2.imag();
              } else {
                SerialWilkinsonShiftInternal::invoke(sub2x2[0], sub2x2[hs1], sub2x2[hs0], sub2x2[hs], &lambda1,
                                                     &lambda2, &is_complex);

                SerialFrancisInternal::invoke(0, mdiff, mdiff, H + hs * mbeg, hs0, hs1, lambda1, lambda2, is_complex);
                /* */ auto &val1    = *(sub2x2 + hs0);
                /* */ auto &val2    = *(sub2x2 - hs1);
                const auto abs_val1 = ats::abs(val1);
                const auto abs_val2 = ats::abs(val2);

                /// convergence check
                if (abs_val1 < tol) {
                  er[(mend - 1) * ers] = sub2x2[hs];
                  ei[(mend - 1) * eis] = zero;
                  val1                 = zero;
                } else if (abs_val2 < tol) {
                  er[(mend - 2) * ers] = lambda1.real();
                  ei[(mend - 2) * eis] = lambda1.imag();
                  er[(mend - 1) * ers] = lambda2.real();
                  ei[(mend - 1) * eis] = lambda2.imag();

                  val1 = zero;
                  val2 = zero;
                }
              }
            }
#endif

          } else {
            /// all eigenvalues are converged
            converge = true;
          }
          ++iter;
        }
        /// Step 3: record missing real eigenvalues from the diagonals
        if (converge) {
          // record undetected eigenvalues
          for (int i = 0; i < m; ++i)
            if (ats::isNan(er[i * ers])) {
              er[i * ers] = H[i * hs];
              ei[i * eis] = zero;
            }
          r_val = 0;
        } else {
          r_val = -1;
        }
        break;
      }
    }
    return r_val;
  }

  /// complex interface
  template <typename RealType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m,
                                           /* */ RealType *H, const int hs0, const int hs1,
                                           /* */ Kokkos::complex<RealType> *e, const int es,
                                           const int max_iteration = 300, const RealType user_tolerence = RealType(-1),
                                           const bool restart = false) {
    RealType *er     = (RealType *)e;
    RealType *ei     = er + 1;
    const int two_es = 2 * es;
    return invoke(m, H, hs0, hs1, er, two_es, ei, two_es, user_tolerence, restart, max_iteration);
  }
};

}  // namespace KokkosBatched

#endif
