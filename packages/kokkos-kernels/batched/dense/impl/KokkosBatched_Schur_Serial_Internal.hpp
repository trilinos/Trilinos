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
#ifndef __KOKKOSBATCHED_SCHUR_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_SCHUR_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_WilkinsonShift_Serial_Internal.hpp"
#include "KokkosBatched_ApplyGivens_Serial_Internal.hpp"
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
struct SerialSchurInternal {
  /// Given a strictly Hessenberg matrix H (m x m), this computes schur
  /// decomposition using the Francis method and stores them into a vector e.
  /// This routine does not scale nor balance the matrix for the numerical
  /// stability.
  ///    H = Z T Z^H and T = Z^H H Z
  /// Parameters:
  ///   [in]m
  ///     A dimension of the square matrix H.
  ///   [in/out]H, [in]hs0, [in]hs1
  ///     Real Hessenberg matrix H(m x m) with strides hs0 and hs1.
  ///     Entering the routine, H is assumed to have a upper Hessenberg form,
  ///     where all subdiagonals are zero. The matrix is overwritten as a upper
  ///     triangular T on exit.
  ///   [in/out]Z, [in]zs0, [in]zs1
  ///     Unitary matrix resulting from Schur decomposition. With a restarting
  ///     option, the matrix may contain previous partial computation results.
  ///   [in/out]w, [in]wlen
  ///     Contiguous workspace of which size is wlen. When restart is true, this
  ///     workspace is not corrupted after the previous iteration. Temporarily,
  ///     it stores subdiag values and given rotations. wlen should be at least
  ///     3*m.
  ///   [in]restart(false)
  ///     With a restart option, the routine assume that the matrix H and the
  ///     vector e contain the partial results from the previous run. When m = 1
  ///     or 2, this option won't work as the routine always computes the all
  ///     eigenvalues.
  ///   [in]user_max_iteration(300)
  ///     Unlike LAPACK which uses various methods for different types of
  ///     matrices, this routine uses the Francis method only. A user can set
  ///     the maximum number of iterations. When it reaches the maximum
  ///     iteration counts without converging all eigenvalues, the routine
  ///     returns -1.
  template <typename RealType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m,
                                           /* */ RealType *H, const int hs0, const int hs1,
                                           /* */ RealType *Z, const int zs0, const int zs1,
                                           /* */ RealType *w, const int wlen, const bool restart = false,
                                           const int user_max_iteration = -1) {
    typedef RealType real_type;
    typedef Kokkos::ArithTraits<real_type> ats;
    const real_type /* one(1), */ zero(0), tol = 1e2 * ats::epsilon();
    const int max_iteration = user_max_iteration < 0 ? 300 : user_max_iteration;
    if (wlen < m * 5) Kokkos::abort("Error: provided workspace is smaller than 3*m");

    int r_val = 0;
    if (restart) {
      if (m <= 2) Kokkos::abort("Error: restart option cannot be used for m=1 or m=2");
    } else {
      /// do not touch input
      /// SerialSetIdentityInternal::invoke(m, Z, zs0, zs1);
    }

    // workspaces
    real_type *subdiags                    = w;
    Kokkos::pair<real_type, real_type> *Gs = (Kokkos::pair<real_type, real_type> *)(w + m);
    if (!restart) {
      /// initialize workspace and Gs
      for (int i = 0; i < m; ++i) subdiags[i] = zero;
    }

    const int hs = hs0 + hs1;  /// diagonal stride
    switch (m) {
      case 0: { /* do nothing */ break;
      }
      case 1: { /* do nothing */ break;
      }
      case 2: {
        /// compute eigenvalues from the characteristic determinant equation
        bool is_complex;
        Kokkos::complex<real_type> lambda1, lambda2;
        Kokkos::pair<real_type, real_type> G;
        SerialSchur2x2Internal::invoke(H, H + hs1, H + hs0, H + hs, &G, &lambda1, &lambda2, &is_complex);

        G.second = -G.second;  // transpose
        SerialApplyRightGivensInternal::invoke(G, 2, Z, zs0, Z + zs1, zs0);
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
          const int mend                    = cnt;
          const int mdiff                   = mend - mbeg;
          const int mend_minus_two_mult_hs0 = (mend - 2) * hs0;

          /// Step 2: if there exist non-converged eigen values
          if (1 < mdiff) {
#if 0  /// implicit QR with shift for testing 
            {
              /// Rayleigh quotient shift 
              const real_type shift = *(H+(mend-1)*hs); 
              SerialHessenbergQR_WithShiftInternal::invoke(mbeg, mend, m, 
                                                           H, hs0, hs1,
                                                           shift,
                                                           Gs, true);

              for (int i=mbeg;i<(mend-1);++i) {
                const Kokkos::pair<real_type,real_type> G(Gs[i].first, -Gs[i].second);
                SerialApplyRightGivensInternal::invoke(G, m,
                                                       Z+i*zs1,     zs0,
                                                       Z+i*zs1+zs1, zs0);                    
              }
            }
#endif

#if 1
            {
              /// find a complex eigen pair
              Kokkos::complex<real_type> lambda1, lambda2;
              bool is_complex;
              real_type *sub2x2 = H + (mend - 2) * hs;
              if (2 == mdiff) {
                Kokkos::pair<real_type, real_type> G;
                SerialSchur2x2Internal::invoke(sub2x2, sub2x2 + hs1, sub2x2 + hs0, sub2x2 + hs, &G, &lambda1, &lambda2,
                                               &is_complex);
                subdiags[mend - 1] = sub2x2[hs0];

                /// apply G' from left
                G.second = -G.second;
                SerialApplyLeftGivensInternal::invoke(G, m - mend, sub2x2 + 2 * hs1, hs1, sub2x2 + hs0 + 2 * hs1, hs1);

                /// apply (G')' from right
                SerialApplyRightGivensInternal::invoke(G, mend - 2, sub2x2 - mend_minus_two_mult_hs0, hs0,
                                                       sub2x2 + hs1 - mend_minus_two_mult_hs0, hs0);
                sub2x2[hs0] = zero;

                /// apply (G')' from right to compute Z
                SerialApplyRightGivensInternal::invoke(G, m, Z + (mend - 2) * zs1, zs0, Z + (mend - 1) * zs1, zs0);

              } else {
                SerialWilkinsonShiftInternal::invoke(sub2x2[0], sub2x2[hs1], sub2x2[hs0], sub2x2[hs], &lambda1,
                                                     &lambda2, &is_complex);

                SerialFrancisInternal::invoke(mbeg, mend, m, H, hs0, hs1, lambda1, lambda2, is_complex, Gs, true);
                /* */ auto &val1    = *(sub2x2 + hs0);
                /* */ auto &val2    = *(sub2x2 - hs1);
                const auto abs_val1 = ats::abs(val1);
                const auto abs_val2 = ats::abs(val2);

                for (int i = mbeg; i < (mend - 1); ++i) {
                  const Kokkos::pair<real_type, real_type> G0(Gs[2 * i].first, -Gs[2 * i].second);
                  const Kokkos::pair<real_type, real_type> G1(Gs[2 * i + 1].first, -Gs[2 * i + 1].second);
                  SerialApplyRightGivensInternal::invoke(G0, m, Z + i * zs1, zs0, Z + i * zs1 + 1 * zs1, zs0);
                  SerialApplyRightGivensInternal::invoke(G1, m, Z + i * zs1, zs0, Z + i * zs1 + 2 * zs1, zs0);
                }

                /// convergence check
                if (abs_val1 < tol) {
                  val1 = zero;
                } else if (abs_val2 < tol) {
                  /// preserve the standard schur form
                  Kokkos::pair<real_type, real_type> G;
                  SerialSchur2x2Internal::invoke(sub2x2, sub2x2 + hs1, sub2x2 + hs0, sub2x2 + hs, &G, &lambda1,
                                                 &lambda2, &is_complex);
                  subdiags[mend - 1] = val1;

                  /// apply G' from left
                  G.second = -G.second;
                  SerialApplyLeftGivensInternal::invoke(G, m - mend, sub2x2 + 2 * hs1, hs1, sub2x2 + hs0 + 2 * hs1,
                                                        hs1);

                  // apply (G')' from right
                  SerialApplyRightGivensInternal::invoke(G, mend - 2, sub2x2 - mend_minus_two_mult_hs0, hs0,
                                                         sub2x2 + hs1 - mend_minus_two_mult_hs0, hs0);
                  val1 = zero;
                  val2 = zero;

                  // apply (G')' from right
                  SerialApplyRightGivensInternal::invoke(G, m, Z + (mend - 2) * zs1, zs0, Z + (mend - 1) * zs1, zs0);
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
          // recover subdiags
          real_type *Hs = H - hs1;
          for (int i = 1; i < m; ++i) {
            Hs[i * hs] = subdiags[i];
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
};

}  // namespace KokkosBatched

#endif
