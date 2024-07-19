// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Details_Lapack128.hpp"
#ifdef HAVE_TEUCHOSCORE_QUADMATH
#  include "Teuchos_BLAS.hpp"
#endif // HAVE_TEUCHOSCORE_QUADMATH


#ifdef HAVE_TEUCHOSCORE_QUADMATH
namespace Teuchos {
namespace Details {

void
Lapack128::
GETRF (const int M, const int N, __float128 A[],
       const int LDA, int IPIV[], int* INFO) const
{
  //std::cerr << "GETRF: N = " << N << std::endl;

  Teuchos::BLAS<int, __float128> blas;

  // NOTE (mfh 05 Sep 2015) This is a direct translation of LAPACK
  // 3.5.0's DGETF2 routine.  LAPACK is under a BSD license.

  *INFO = 0;
  if (M < 0) {
    *INFO = -1;
  } else if (N < 0) {
    *INFO = -2;
  } else if (LDA < std::max (1, M)) {
    *INFO = -4;
  }
  if (*INFO != 0) {
    return;
  }

  // Quick return if possible
  if (M == 0 || N == 0) {
    return;
  }

  // Compute machine safe minimum sfmin (such that 1/sfmin does
  // not overflow).  LAPACK 3.1 just returns for this the smallest
  // normalized number.
  const __float128 sfmin = FLT128_MIN;
  const __float128 zero = 0.0;
  const __float128 one = 1.0;

  const int j_upperBound = std::min (M, N);
  for (int j = 1; j <= j_upperBound; ++j) {
    //std::cerr << "  j = " << j << std::endl;

    // Find pivot and test for singularity.
    __float128* const A_jj = A + (j-1)*LDA + (j-1);

    //std::cerr << "  CALLING IAMAX" << std::endl;
    const int jp = (j - 1) + blas.IAMAX (M - j + 1, A_jj, 1);
    IPIV[j - 1] = jp;

    const __float128* A_jp_j = A + jp + LDA*j;
    if (*A_jp_j != zero) {
      // Apply the interchange to columns 1:N.
      __float128* const A_j1 = A + (j - 1);
      __float128* const A_jp_1 = A + (jp - 1);

      if (jp != j) {
        blas.SWAP (N, A_j1, LDA, A_jp_1, LDA);
      }

      // Compute elements J+1:M of J-th column.
      if (j < M) {
            __float128* const A_j1_j = A + j + (j-1)*LDA;

        if (fabsq (*A_jj) >= sfmin) {
          blas.SCAL (M-j, one / *A_jj, A_j1_j, 1);
        } else {
          for (int i = 1; i <= M-j; ++i) {
            __float128* const A_jpi_j = A + (j+i-1) + (j-1)*LDA;
            *A_jpi_j /= *A_jj;
          }
        }
      }
    } else if (*INFO == 0) {
      *INFO = j;
    }

    if (j < std::min (M, N)) {
      //std::cerr << "  UPDATE TRAILING SUBMATRIX" << std::endl;

      // Update trailing submatrix.
      const __float128* A_j1_j = A + j + (j-1)*LDA;
      const __float128* A_j_j1 = A + (j-1) + j*LDA;
      __float128* A_j1_j1 = A + j + j*LDA;
      blas.GER (M-j, N-j, -one, A_j1_j, 1, A_j_j1, LDA, A_j1_j1, LDA);
    }
  }
}

void
Lapack128::
LASWP (const int N, __float128 A[], const int LDA, const int K1,
       const int K2, const int IPIV[], const int INCX) const
{
  int i, i1, i2, inc, ip, ix, ix0, j, k, n32;
  __float128 temp;

  // Interchange row I with row IPIV(I) for each of rows K1 through K2.

  if (INCX > 0) {
    ix0 = K1;
    i1 = K1;
    i2 = K2;
    inc = 1;
  } else if (INCX < 0) {
    ix0 = 1 + (1 - K2)*INCX;
    i1 = K2;
    i2 = K1;
    inc = -1;
  } else { // INCX == 0
    return;
  }

  // The LAPACK 3.5.0 source code does 32 entries at a time,
  // presumably for better vectorization or cache line usage.
  n32 = (N / 32) * 32;

  if (n32 != 0) {
    for (j = 1; j <= n32; j += 32) {
      ix = ix0;
      // C and C++ lack Fortran's convenient range specifier,
      // which can iterate over a range in either direction
      // without particular fuss about the end condition.
      for (i = i1; (inc > 0) ? (i <= i2) : (i >= i2); i += inc) {
        ip = IPIV[ix-1];
        if (ip != i) {
          for (k = j; k <= j+31; ++k) {
            temp = A[(i-1) + (k-1)*LDA]; //  temp = a( i, k )
            A[(i-1) + (k-1)*LDA] = A[(ip-1) + (k-1)*LDA]; // a( i, k ) = a( ip, k )
            A[(ip-1) + (k-1)*LDA] = temp; // a( ip, k ) = temp
          }
        }
        ix = ix + INCX;
      }
    }
  }

  if (n32 != N) {
    n32 = n32 + 1;
    ix = ix0;
    // C and C++ lack Fortran's convenient range specifier,
    // which can iterate over a range in either direction
    // without particular fuss about the end condition.
    for (i = i1; (inc > 0) ? (i <= i2) : (i >= i2); i += inc) {
      ip = IPIV[ix-1];
      if (ip != i) {
        for (k = n32; k <= N; ++k) {
          temp = A[(i-1) + (k-1)*LDA]; //  temp = a( i, k )
          A[(i-1) + (k-1)*LDA] = A[(ip-1) + (k-1)*LDA]; // a( i, k ) = a( ip, k )
          A[(ip-1) + (k-1)*LDA] = temp; // a( ip, k ) = temp
        }
      }
      ix = ix + INCX;
    }
  }
}

void
Lapack128::
GETRI (const int /* N */, __float128 /* A */ [], const int /* LDA */,
       int /* IPIV */ [], __float128 /* WORK */ [], const int /* LWORK */,
       int* /* INFO */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, __float128>::GETRI: Not implemented yet.");
}


void
Lapack128::
GETRS (const char TRANS, const int N, const int NRHS,
       const __float128 A[], const int LDA, const int IPIV[],
       __float128 B[], const int LDB, int* INFO) const
{
  //std::cerr << "GETRS: N = " << N << std::endl;

  Teuchos::BLAS<int, __float128> blas;

  // NOTE (mfh 05 Sep 2015) This is a direct translation of LAPACK
  // 3.5.0's DGETRS routine.  LAPACK is under a BSD license.

  *INFO = 0;
  const bool notran = (TRANS == 'N' || TRANS == 'n');
  if (! notran
      && ! (TRANS == 'T' || TRANS == 't')
      && ! (TRANS == 'C' || TRANS == 'c')) {
    *INFO = -1; // invalid TRANS argument
  }
  else if (N < 0) {
    *INFO = -2; // invalid N (negative)
  }
  else if (NRHS < 0) {
    *INFO = -3; // invalid NRHS (negative)
  }
  else if (LDA < std::max (1, N)) {
    *INFO = -5; // invalid LDA (too small)
  }
  else if (LDB < std::max (1, N)) {
    *INFO = -8; // invalid LDB (too small)
  }
  if (*INFO != 0) {
    return;
  }

  const __float128 one = 1.0;

  using Teuchos::LEFT_SIDE;
  using Teuchos::LOWER_TRI;
  using Teuchos::UPPER_TRI;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;

  if (notran) { // No transpose; solve AX=B
    // Apply row interchanges to the right-hand sides.
    //std::cerr << "CALLING LASWP" << std::endl;
    LASWP (NRHS, B, LDB, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    //std::cerr << "CALLING TRSM (1)" << std::endl;
    blas.TRSM (LEFT_SIDE, LOWER_TRI, NO_TRANS, UNIT_DIAG, N, NRHS,
               one, A, LDA, B, LDB);
    // Solve U*X = B, overwriting B with X.
    //std::cerr << "CALLING TRSM (2)" << std::endl;
    blas.TRSM (LEFT_SIDE, UPPER_TRI, NO_TRANS, NON_UNIT_DIAG, N, NRHS,
               one, A, LDA, B, LDB);
  }
  else { // Transpose or conjugate transpose: solve A^{T,H}X = B.
    const Teuchos::ETransp transposeMode = (TRANS == 'T' || TRANS == 't') ?
      TRANS : CONJ_TRANS;

    // Solve U^{T,H}*X = B, overwriting B with X.
    //std::cerr << "CALLING TRSM (1)" << std::endl;
    blas.TRSM (LEFT_SIDE, UPPER_TRI, transposeMode, NON_UNIT_DIAG, N, NRHS,
               one, A, LDA, B, LDB);
    // Solve L^{T,H}*X = B, overwriting B with X.
    //std::cerr << "CALLING TRSM (2)" << std::endl;
    blas.TRSM (LEFT_SIDE, LOWER_TRI, transposeMode, UNIT_DIAG, N, NRHS,
               one, A, LDA, B, LDB);
    //std::cerr << "CALLING LASWP" << std::endl;
    // Apply row interchanges to the solution vectors.
    LASWP (NRHS, B, LDB, 1, N, IPIV, -1);
  }

  //std::cerr << "DONE WITH GETRS" << std::endl;
}

__float128
Lapack128::
LAPY2 (const __float128& x, const __float128& y) const
{
  const __float128 xabs = fabsq (x);
  const __float128 yabs = fabsq (y);
  const __float128 w = fmaxq (xabs, yabs);
  const __float128 z = fminq (xabs, yabs);

  if (z == 0.0) {
    return w;
  } else {
    const __float128 one = 1.0;
    const __float128 z_div_w = z / w;
    return w * sqrtq (one + z_div_w * z_div_w);
  }
}

void
Lapack128::
ORM2R (const char side, const char trans,
       const int m, const int n, const int k,
       const __float128 A[], const int lda,
       const __float128* const tau,
       __float128 C[], const int ldc,
       __float128 work[], int* const info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}

namespace { // (anonymous)

  int
  ILADLC (const int m, const int n, const __float128 A[], const int lda)
  {
    const __float128 zero = 0.0;

    // Quick test for the common case where one corner is non-zero.
    if (n == 0) {
      return n;
    } else if (A[0 + (n-1)*lda] != zero || A[(m-1) + (n-1)*lda] != zero) {
      return n;
    } else {
      // Now scan each column from the end, returning with the first non-zero.
      for (int j = n; j > 0; --j) {
        for (int i = 1; i <= m; ++i) {
          if (A[(i-1) + (j-1)*lda] != zero) {
            return j;
          }
        }
      }
      return 0;
    }
  }

  int
  ILADLR (const int m, const int n, const __float128 A[], const int lda)
  {
    const __float128 zero = 0.0;

    // Quick test for the common case where one corner is non-zero.
    if (m == 0) {
      return m;
    } else if (A[(m-1) + 0*lda] != zero || A[(m-1) + (n-1)*lda] != zero) {
      return m;
    } else {
      // Scan up each column tracking the last zero row seen.
      int lastZeroRow = 0;
      for (int j = 1; j <= n; ++j) {
        int i = m;
        while (A[(std::max (i, 1) - 1) + (j - 1)*lda] == zero && i >= 1) {
          i--;
        }
        lastZeroRow = std::max (lastZeroRow, i);
      }
      return lastZeroRow;
    }
  }
} // namespace (anonymous)

void
Lapack128::
LARF (const char side,
      const int m,
      const int n,
      const __float128 v[],
      const int incv,
      const __float128 tau,
      __float128 C[],
      const int ldc,
      __float128 work[]) const
{
  const __float128 zero = 0.0;
  const __float128 one = 1.0;
  Teuchos::BLAS<int, __float128> blas;
  const bool applyLeft = (side == 'L');
  int lastv = 0;
  int lastc = 0;
  int i = 0;

  if (tau != zero) {
    // Set up variables for scanning V.  LASTV begins pointing to the end of V.
    if (applyLeft) {
      lastv = m;
    } else {
      lastv = n;
    }
    if (incv > 0) {
      i = 1 + (lastv - 1) * incv;
    } else {
      i = 1;
    }
    // Look for the last non-zero row in V.
    while (lastv > 0 && v[i-1] == zero) {
      lastv = lastv - 1;
      i = i - incv;
    }
    if (applyLeft) {
      // Scan for the last non-zero column in C(1:lastv,:).
      lastc = ILADLC (lastv, n, C, ldc);
    } else {
      // Scan for the last non-zero row in C(:,1:lastv).
      lastc = ILADLR (m, lastv, C, ldc);
    }
  }

  // Note that lastc == 0 renders the BLAS operations null; no special
  // case is needed at this level.
  if (applyLeft) {
    // Form  H * C
    if (lastv > 0) {
      // w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
      blas.GEMV (Teuchos::TRANS, lastv, lastc, one, C, ldc, v, incv,
                 zero, work, 1);
      // C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
      blas.GER (lastv, lastc, -tau, v, incv, work, 1, C, ldc);
    }
  }
  else {
    // Form  C * H
    if (lastv > 0) {
      // w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
      blas.GEMV (Teuchos::NO_TRANS, lastc, lastv, one, C, ldc,
                 v, incv, zero, work, 1);
      // C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
      blas.GER (lastc, lastv, -tau, work, 1, v, incv, C, ldc);
    }
  }
}


void
Lapack128::
LARFG (const int N, __float128* const ALPHA,
       __float128 X[], const int INCX, __float128* const TAU) const
{
  // This is actually LARFGP.

  const __float128 zero = 0.0;
  const __float128 one = 1.0;
  const __float128 two = 2.0;
  Teuchos::BLAS<int, __float128> blas;

  if (N <= 0) {
    *TAU = zero;
    return;
  }
  __float128 xnorm = blas.NRM2 (N-1, X, INCX);

  if (xnorm == zero) {
    // H  =  [+/-1, 0; I], sign chosen so *ALPHA >= 0
    if (*ALPHA >= zero) {
      // When TAU.eq.ZERO, the vector is special-cased to be all zeros
      // in the application routines.  We do not need to clear it.
      *TAU = zero;
    } else {
      // However, the application routines rely on explicit
      // zero checks when TAU.ne.ZERO, and we must clear X.
      *TAU = two;
      for (int j = 0; j < N; ++j) {
        X[j * INCX] = 0.0;
      }
      *ALPHA = -*ALPHA;
    }
  } else { // general case (norm of x is nonzero)
    // This implements Fortran's two-argument SIGN intrinsic.
    __float128 beta = copysignq (LAPY2 (*ALPHA, xnorm), *ALPHA);
    const __float128 smlnum = FLT128_MIN / FLT128_EPSILON;
    int knt = 0;

    if (fabsq (beta) < smlnum) {
      // XNORM, BETA may be inaccurate; scale X and recompute them

      __float128 bignum = one / smlnum;
      do {
        knt = knt + 1;
        blas.SCAL (N-1, bignum, X, INCX);
        beta = beta*bignum;
        *ALPHA = *ALPHA*bignum;
      } while (fabsq(beta) < smlnum);

      // New BETA is at most 1, at least SMLNUM
      xnorm = blas.NRM2 (N-1, X, INCX);
      beta = copysignq (LAPY2 (*ALPHA, xnorm), *ALPHA);
    }

    __float128 savealpha = *ALPHA;
    *ALPHA = *ALPHA + beta;
    if (beta < zero) {
      beta = -beta;
      *TAU = -*ALPHA / beta;
    } else {
      *ALPHA = xnorm * (xnorm / *ALPHA);
      *TAU = *ALPHA / beta;
      *ALPHA = -*ALPHA;
    }

    if (fabsq (*TAU) <= smlnum) {
      // In the case where the computed TAU ends up being a
      // denormalized number, it loses relative accuracy. This is a
      // BIG problem. Solution: flush TAU to ZERO. This explains the
      // next IF statement.
      //
      // (Bug report provided by Pat Quillen from MathWorks on Jul 29,
      // 2009.)  (Thanks Pat. Thanks MathWorks.)

      if (savealpha >= zero) {
        *TAU = zero;
      } else {
        *TAU = two;
        for (int j = 0; j < N; ++j) {
          X[j*INCX] = 0.0;
        }
        beta = -savealpha;
      }
    }
    else { // this is the general case
      blas.SCAL (N-1, one / *ALPHA, X, INCX);
    }
    // If BETA is subnormal, it may lose relative accuracy
    for (int j = 1; j <= knt; ++j) {
      beta = beta*smlnum;
    }
    *ALPHA = beta;
  }
}

void
Lapack128::
GEQR2 (const int /* M */,
       const int /* N */,
       __float128 /* A */ [],
       const int /* LDA */,
       __float128 /* TAU */ [],
       __float128 /* WORK */ [],
       int* const /* INFO */ ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, __float128>::GEQR2: Not implemented yet.");
}

void
Lapack128::
GEQRF (const int M,
       const int N,
       __float128 A[],
       const int LDA,
       __float128 TAU[],
       __float128 WORK[],
       const int LWORK,
       int* const INFO) const
{
  // mfh 17 Sep 2015: We don't implement a BLAS 3 QR factorization for
  // __float128.  Instead, we call the BLAS 2 QR factorization GEQR2,
  // which has a fixed minimum WORK array length of N.  Thus, we have
  // to roll our own LWORK query here.

  if (LWORK == -1) {
    WORK[0] = static_cast<__float128> (N);
  }
  else {
    GEQR2 (M, N, A, LDA, TAU, WORK, INFO);
  }
}

void
Lapack128::
ORGQR (const int /* M */,
       const int /* N */,
       const int /* K */,
       __float128 /* A */ [],
       const int /* LDA */,
       const __float128 /* TAU */ [],
       __float128 /* WORK */ [],
       const int /* LWORK */,
       int* const /* INFO */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, __float128>::GEQR2: Not implemented yet.");
}

void
Lapack128::
UNGQR (const int /* M */,
       const int /* N */,
       const int /* K */,
       __float128 /* A */ [],
       const int /* LDA */,
       const __float128 /* TAU */ [],
       __float128 /* WORK */ [],
       const int /* LWORK */,
       int* const /* INFO */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, __float128>::GEQR2: Not implemented yet.");
}

void
Lapack128::
LASCL (const char TYPE,
       const int kl,
       const int ku,
       const __float128 cfrom,
       const __float128 cto,
       const int m,
       const int n,
       __float128* A,
       const int lda,
       int* info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, __float128>::LASCL: Not implemented yet.");
}

void
Lapack128::
GBTRF (const int m,
       const int n,
       const int kl,
       const int ku,
       __float128* A,
       const int lda,
       int* IPIV,
       int* info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, __float128>::GBTRF: Not implemented yet.");
}

void
Lapack128::
GBTRS (const char TRANS,
       const int n,
       const int kl,
       const int ku,
       const int nrhs,
       const __float128* A,
       const int lda,
       const int* IPIV,
       __float128* B,
       const int ldb,
       int* info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, __float128>::GBTRS: Not implemented yet.");
}

} // namespace Details
} // namespace Teuchos
#endif // HAVE_TEUCHOSCORE_QUADMATH
