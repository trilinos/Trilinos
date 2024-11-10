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
#ifndef __KOKKOSBATCHED_SVD_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_SVD_SERIAL_INTERNAL_HPP__

/// \author Brian Kelley (bmkelle@sandia.gov)

#include "Kokkos_MathematicalFunctions.hpp"
#include "KokkosBatched_SetIdentity_Internal.hpp"
#include "KokkosBatched_Givens_Serial_Internal.hpp"
#include "KokkosBatched_ApplyGivens_Serial_Internal.hpp"
#include "KokkosBatched_Householder_Serial_Internal.hpp"
#include "KokkosBatched_ApplyHouseholder_Serial_Internal.hpp"

// Use this macro to handle raw pointer/stride based 2D indexing in this file
// (just for readability) Requires that for pointer X, the corresponding row/col
// strides are named Xs0 and Xs1.
#define SVDIND(arr, i, j) arr[(i)*arr##s0 + (j)*arr##s1]
#define SVDSWAP(a, b) \
  {                   \
    auto tmp = a;     \
    a        = b;     \
    b        = tmp;   \
  }

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

struct SerialSVDInternal {
  // Find the two eigenvalues of [a11 a21 ; a21 a22] by solving the
  // characteristic quadratic. Since matrix is symmetric these will be real.
  // NOTE: this is essentially the Wilkinson shift routine already in Batched,
  // however this is simpler because it exploits the symmetric structure, and
  // the realness of the eigenvalues.
  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static void symEigen2x2(value_type a11, value_type a21, value_type a22, value_type& e1,
                                                 value_type& e2) {
    value_type a       = Kokkos::ArithTraits<value_type>::one();
    value_type b       = -a11 - a22;
    value_type c       = a11 * a22 - a21 * a21;
    value_type sqrtDet = Kokkos::sqrt(b * b - 4 * a * c);
    e1                 = (-b + sqrtDet) / (2 * a);
    e2                 = (-b - sqrtDet) / (2 * a);
  }

  // B is a square submatrix on the diagonal.
  // Usub is a subset of columns of U
  // Vtsub is a subset of rows of Vt
  //
  // B22 is nsub * nsub, Usub is m * nsub, and Vtsub is nsub * n
  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static void svdStep(value_type* B, value_type* U, value_type* Vt, int um, int vn, int n,
                                             int Bs0, int Bs1, int Us0, int Us1, int Vts0, int Vts1) {
    using KAT = Kokkos::ArithTraits<value_type>;
    // Compute the eigenvalues of trailing 2x2
    value_type dn     = SVDIND(B, n - 1, n - 1);
    value_type dm     = SVDIND(B, n - 2, n - 2);
    value_type fm     = SVDIND(B, n - 2, n - 1);
    value_type fmm1   = (n > 2) ? SVDIND(B, n - 3, n - 2) : KAT::zero();
    value_type target = dn * dn + fm * fm;
    value_type e1, e2, mu;
    symEigen2x2(dm * dm + fmm1 * fmm1, dm * fm, target, e1, e2);
    // the shift is the eigenvalue closer to the last diagonal entry of B^T*B
    if (Kokkos::abs(e1 - target) < Kokkos::abs(e2 - target))
      mu = e1;
    else
      mu = e2;
    value_type y = SVDIND(B, 0, 0) * SVDIND(B, 0, 0) - mu;
    value_type z = SVDIND(B, 0, 0) * SVDIND(B, 0, 1);
    for (int k = 0; k < n - 1; k++) {
      // Use Givens to zero out z in [y; z]
      Kokkos::pair<value_type, value_type> G;
      value_type discard;  // Don't actually write [alpha; 0] anywhere
      KokkosBatched::SerialGivensInternal::invoke<value_type>(y, z, &G, &discard);
      // apply the Givens transformation to B on the right, to columns k,k+1
      // B := BG(k, k+1, theta)
      int minrow = KOKKOSKERNELS_MACRO_MAX(0, k - 1);
      int maxrow = KOKKOSKERNELS_MACRO_MIN(n, k + 2);
      KokkosBatched::SerialApplyRightGivensInternal::invoke<value_type>(G, maxrow - minrow, &SVDIND(B, minrow, k + 1),
                                                                        Bs0, &SVDIND(B, minrow, k), Bs0);
      if (Vt) {
        KokkosBatched::SerialApplyLeftGivensInternal::invoke<value_type>(G, vn, &SVDIND(Vt, k + 1, 0), Vts1,
                                                                         &SVDIND(Vt, k, 0), Vts1);
      }
      y = SVDIND(B, k, k);
      z = SVDIND(B, k + 1, k);
      KokkosBatched::SerialGivensInternal::invoke<value_type>(y, z, &G, &SVDIND(B, k, k));
      SVDIND(B, k + 1, k) = KAT::zero();
      int mincol          = k + 1;
      int maxcol          = KOKKOSKERNELS_MACRO_MIN(n, k + 3);
      // apply Givens transformation to B on the left, to rows k, k + 1
      // B := G(k, k+1, theta)^T * B
      KokkosBatched::SerialApplyLeftGivensInternal::invoke<value_type>(G, maxcol - mincol, &SVDIND(B, k + 1, mincol),
                                                                       Bs1, &SVDIND(B, k, mincol), Bs1);
      if (U) {
        KokkosBatched::SerialApplyRightGivensInternal::invoke<value_type>(G, um, &SVDIND(U, 0, k + 1), Us0,
                                                                          &SVDIND(U, 0, k), Us0);
      }
      if (k < n - 2) {
        y = SVDIND(B, k, k + 1);
        z = SVDIND(B, k, k + 2);
      }
    }
  }

  // Deal with B(i, i) = 0, by chasing superdiagonal nonzero across row i.
  // Assumes i is not the last row.
  // U is m*m, B is n*n
  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static void svdZeroRow(int i, value_type* B, int n, int Bs0, int Bs1, value_type* U, int Um,
                                                int Us0, int Us1) {
    Kokkos::pair<value_type, value_type> G;
    for (int j = i + 1; j < n; j++) {
      // Zero out B(i, j) against diagonal j, introducing nonzero in B(i, j + 1)
      KokkosBatched::SerialGivensInternal::invoke<value_type>(SVDIND(B, j, j), SVDIND(B, i, j), &G, &SVDIND(B, j, j));
      SVDIND(B, i, j) = Kokkos::ArithTraits<value_type>::zero();
      // Now, only need to apply givens to a single column (if not already at
      // the end), introducing the next nonzero
      if (j < n - 1) {
        KokkosBatched::SerialApplyLeftGivensInternal::invoke<value_type>(G, 1, &SVDIND(B, i, j + 1), Bs1,
                                                                         &SVDIND(B, j, j + 1), Bs1);
      }
      if (U) {
        KokkosBatched::SerialApplyRightGivensInternal::invoke<value_type>(G, Um, &SVDIND(U, 0, i), Us0,
                                                                          &SVDIND(U, 0, j), Us0);
      }
    }
  }

  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static void svdZeroLastColumn(value_type* B, int n, int Bs0, int Bs1, int vn, value_type* Vt,
                                                       int Vts0, int Vts1) {
    // Deal with B(n-1, n-1) = 0, by chasing the superdiagonal nonzero up the last column.
    Kokkos::pair<value_type, value_type> G;
    for (int j = n - 2; j >= 0; j--) {
      KokkosBatched::SerialGivensInternal::invoke<value_type>(SVDIND(B, j, j), SVDIND(B, j, n - 1), &G,
                                                              &SVDIND(B, j, j));
      SVDIND(B, j, n - 1) = Kokkos::ArithTraits<value_type>::zero();
      if (j != 0) {
        KokkosBatched::SerialApplyRightGivensInternal::invoke<value_type>(G, 1, &SVDIND(B, j - 1, n - 1), Bs0,
                                                                          &SVDIND(B, j - 1, j), Bs0);
      }
      if (Vt) {
        KokkosBatched::SerialApplyLeftGivensInternal::invoke<value_type>(G, vn, &SVDIND(Vt, n - 1, 0), Vts1,
                                                                         &SVDIND(Vt, j, 0), Vts1);
      }
    }
  }

  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static void bidiagonalize(int m, int n, value_type* A, int As0, int As1, value_type* U,
                                                   int Us0, int Us1, value_type* Vt, int Vts0, int Vts1,
                                                   value_type* work) {
    using KAT = Kokkos::ArithTraits<value_type>;
    value_type tau;
    for (int i = 0; i < n; i++) {
      // Eliminating column i of A below the diagonal
      KokkosBatched::SerialLeftHouseholderInternal::invoke<value_type>(m - i - 1, &SVDIND(A, i, i),
                                                                       &SVDIND(A, i + 1, i), As0, &tau);
      if (n - i > 1) {
        KokkosBatched::SerialApplyLeftHouseholderInternal::invoke<value_type>(
            m - i - 1, n - i - 1, &tau, &SVDIND(A, i + 1, i), As0, &SVDIND(A, i, i + 1), As1, &SVDIND(A, i + 1, i + 1),
            As0, As1, work);
      }
      if (U) {
        KokkosBatched::SerialApplyRightHouseholderInternal::invoke<value_type>(
            m, m - i - 1, &tau, &SVDIND(A, i + 1, i), As0, &SVDIND(U, 0, i), Us0, &SVDIND(U, 0, i + 1), Us0, Us1, work);
      }
      // Zero out A subdiag explicitly (NOTE: may not be necessary...)
      for (int j = i + 1; j < m; j++) {
        SVDIND(A, j, i) = KAT::zero();
      }
      if (i < n - 2) {
        // Eliminating row i of A to the right of the 1st superdiagonal
        KokkosBatched::SerialLeftHouseholderInternal::invoke<value_type>(n - i - 2, &SVDIND(A, i, i + 1),
                                                                         &SVDIND(A, i, i + 2), As1, &tau);
        if (m - i > 1) {
          KokkosBatched::SerialApplyRightHouseholderInternal::invoke<value_type>(
              m - i - 1, n - i - 2, &tau, &SVDIND(A, i, i + 2), As1, &SVDIND(A, i + 1, i + 1), As0,
              &SVDIND(A, i + 1, i + 2), As0, As1, work);
        }
        if (Vt) {
          KokkosBatched::SerialApplyLeftHouseholderInternal::invoke<value_type>(
              n - i - 2, n, &tau, &SVDIND(A, i, i + 2), As1, &SVDIND(Vt, i + 1, 0), Vts1, &SVDIND(Vt, i + 2, 0), Vts0,
              Vts1, work);
        }
        // Zero out A superdiag row explicitly
        for (int j = i + 2; j < n; j++) {
          SVDIND(A, i, j) = KAT::zero();
        }
      }
    }
  }

  // Compute the SVD of a bidiagonal matrix B. Apply inverse transformations to
  // U and Vt to maintain the product U*B*Vt. At the end, the singular values
  // are copied to sigma.
  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static void bidiSVD(int m, int n, value_type* B, int Bs0, int Bs1, value_type* U, int Us0,
                                             int Us1, value_type* Vt, int Vts0, int Vts1, value_type* sigma, int ss,
                                             const value_type& tol) {
    using KAT            = Kokkos::ArithTraits<value_type>;
    const value_type eps = Kokkos::ArithTraits<value_type>::epsilon();
    int p                = 0;
    int q                = 0;
    while (true) {
      // Zero out tiny superdiagonal entries
      for (int i = 0; i < n - 1; i++) {
        if (Kokkos::abs(SVDIND(B, i, i + 1)) <
                eps * (Kokkos::abs(SVDIND(B, i, i)) + Kokkos::abs(SVDIND(B, i + 1, i + 1))) ||
            Kokkos::abs(SVDIND(B, i, i + 1)) < tol) {
          SVDIND(B, i, i + 1) = KAT::zero();
        }
      }
      // Find q: first column from the end with nonzero superdiagonal.
      // If no such columns, will be 0.
      for (q = n - 1; q > 0; q--) {
        if (SVDIND(B, q - 1, q) != KAT::zero()) break;
      }
      if (q == 0) {
        // B is completely diagonal, so it contains singular values and we are
        // done.
        break;
      }
      q++;
      // now, q is the upper (exclusive) bound of submatrix on which to do SVD
      // step. Find min p, so that [p, q) x [p, q) submatrix has all nonzero
      // superdiagonals.
      for (p = q - 1; p > 0; p--) {
        if (SVDIND(B, p - 1, p) == KAT::zero()) break;
      }
      value_type* Bsub  = &SVDIND(B, p, p);
      value_type* Usub  = &SVDIND(U, 0, p);
      value_type* Vtsub = &SVDIND(Vt, p, 0);
      int nsub          = q - p;
      // If there are zero diagonals in this range, eliminate the entire row
      //(effectively decoupling into two subproblems)
      for (int i = q - 1; i >= p; i--) {
        if (SVDIND(B, i, i) == KAT::zero()) {
          if (i == q - 1) {
            // Last diagonal entry being 0 is a special case.
            // Zero out the superdiagonal above it.
            // Deal with B(q-1, q-1) = 0, by chasing the superdiagonal nonzero
            // B(q-2, q-1) up the last column.
            //
            // Once that nonzero reaches B(p, q-1), we are either at the top of B
            // (if p == 0) or the superdiag above B(p, p) is zero.
            // In either case, the chase stops after eliminating B(p, q-1) because no
            // new entry is introduced by the Givens.
            svdZeroLastColumn(Bsub, nsub, Bs0, Bs1, n, Vtsub, Vts0, Vts1);
          } else if (SVDIND(B, i, i + 1) != KAT::zero()) {
            svdZeroRow(i - p, Bsub, nsub, Bs0, Bs1, Usub, m, Us0, Us1);
          }
        }
      }
      // B22 is nsub * nsub, Usub is m * nsub, and Vtsub is nsub * n
      svdStep(Bsub, Usub, Vtsub, m, n, nsub, Bs0, Bs1, Us0, Us1, Vts0, Vts1);
    }
    for (int i = 0; i < n; i++) {
      sigma[i * ss] = SVDIND(B, i, i);
    }
  }

  // Convert SVD into conventional form: singular values positive and in
  // descending order
  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static void postprocessSVD(int m, int n, value_type* U, int Us0, int Us1, value_type* Vt,
                                                    int Vts0, int Vts1, value_type* sigma, int ss) {
    // First step: flip signs on negative singular values
    for (int i = 0; i < n; i++) {
      if (sigma[i * ss] < 0) {
        sigma[i * ss] = -sigma[i * ss];
        if (Vt) {
          for (int j = 0; j < n; j++) SVDIND(Vt, i, j) = -SVDIND(Vt, i, j);
        }
      }
    }
    // Second step: stable selection sort to put singular values in order.
    // Using selection sort because the quadratic part only applies to sigma
    // (O(n^2) total), and it minimizes column swaps in U,V (O(mn) total
    // movement).
    for (int i = 0; i < n - 1; i++) {
      // find the proper singular value to go in position i
      value_type maxval = sigma[i * ss];
      int maxloc        = i;
      for (int j = i + 1; j < n; j++) {
        if (sigma[j * ss] > maxval) {
          maxval = sigma[j * ss];
          maxloc = j;
        }
      }
      // swap singular values and U/V columns i and maxloc (if maxloc is not
      // already in the right place)
      if (i != maxloc) {
        SVDSWAP(sigma[i * ss], sigma[maxloc * ss]);
        if (U) {
          for (int j = 0; j < m; j++) SVDSWAP(SVDIND(U, j, i), SVDIND(U, j, maxloc))
        }
        if (Vt) {
          for (int j = 0; j < n; j++) SVDSWAP(SVDIND(Vt, i, j), SVDIND(Vt, maxloc, j))
        }
      }
    }
  }

  template <typename value_type>
  KOKKOS_INLINE_FUNCTION static int invoke(int m, int n, value_type* A, int As0, int As1, value_type* U, int Us0,
                                           int Us1, value_type* Vt, int Vts0, int Vts1, value_type* sigma, int ss,
                                           value_type* work, value_type tol = Kokkos::ArithTraits<value_type>::zero()) {
    // First, if m < n, need to instead compute (V, s, U^T) = A^T.
    // This just means swapping U & Vt, and implicitly transposing A, U and Vt.
    if (m < n) {
      // Transpose A
      SVDSWAP(m, n);
      SVDSWAP(As0, As1);
      // Transpose and swap U, Vt
      SVDSWAP(U, Vt);
      SVDSWAP(Us0, Vts1);
      SVDSWAP(Us1, Vts0);
    }
    if (U) {
      KokkosBatched::SerialSetIdentityInternal::invoke<value_type>(m, m, U, Us0, Us1);
    }
    if (Vt) {
      KokkosBatched::SerialSetIdentityInternal::invoke<value_type>(n, n, Vt, Vts0, Vts1);
    }
    if (m == 0 || n == 0) {
      // sigma is length 0, so there's nothing left to compute
      return 0;
    }
    bidiagonalize(m, n, A, As0, As1, U, Us0, Us1, Vt, Vts0, Vts1, work);
    bidiSVD(m, n, A, As0, As1, U, Us0, Us1, Vt, Vts0, Vts1, sigma, ss, tol);
    postprocessSVD(m, n, U, Us0, Us1, Vt, Vts0, Vts1, sigma, ss);
    return 0;
  }
};

}  // namespace KokkosBatched

#undef SVDIND
#undef SVDSWAP

#endif
