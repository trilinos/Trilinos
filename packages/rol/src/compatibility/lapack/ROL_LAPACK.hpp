// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_LAPACK_H
#define ROL_LAPACK_H

/** \class ROL::LAPACK
  \brief Provides interface to Lapack
  */

/* A) Define PREFIX and Teuchos_fcd based on platform. */


// Following two macros were stolen from TEUCHOS_LAPACK.cpp
/* for INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
*/
#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, one
#else
#define CHAR_MACRO(char_var) &char_var
#endif

#ifdef CHARPTR_MACRO
#undef CHARPR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHARPTR_MACRO(charptr_var) charptr_var, one
#else
#define CHARPTR_MACRO(charptr_var) charptr_var
#endif

// Following macros were stolen from TEUCHOS_LAPACK_wrappers.hpp
#if defined(INTEL_CXML)
#  define PREFIX __stdcall
#  define Teuchos_fcd const char *, unsigned int
#elif defined(INTEL_MKL)
#  define PREFIX
#  define Teuchos_fcd const char *
#else
#  define PREFIX
#  define Teuchos_fcd const char *
#endif

// in progress - added for EQUED which is a modified char *, not const
#if defined(INTEL_CXML)
#  define PREFIX __stdcall
#  define Teuchos_nonconst_fcd char *, unsigned int // Need to evaluate unsigned int - CXML deprecated
#elif defined(INTEL_MKL)
#  define PREFIX
#  define Teuchos_nonconst_fcd char *
#else
#  define PREFIX
#  define Teuchos_nonconst_fcd char *
#endif

#define DGEEV_F77  F77_BLAS_MANGLE(dgeev, DGEEV)
#define DLATRS_F77  F77_BLAS_MANGLE(dlatrs, DLATRS)
#define DGTTRS_F77  F77_BLAS_MANGLE(dgttrs, DGTTRS)
#define DGTTRF_F77  F77_BLAS_MANGLE(dgttrf, DGTTRF)
#define DTRTRS_F77  F77_BLAS_MANGLE(dtrtrs, DTRTRS)
#define DGETRF_F77  F77_BLAS_MANGLE(dgetrf, DGETRF)
#define DGETRI_F77  F77_BLAS_MANGLE(dgetri, DGETRI)
#define DGETRS_F77  F77_BLAS_MANGLE(dgetrs, DGETRS)
#define DGEQRF_F77  F77_BLAS_MANGLE(dgeqrf, DGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(dorgqr, DORGQR)
#define DPTTRS_F77  F77_BLAS_MANGLE(dpttrs, DPTTRS)
#define DPTTRF_F77  F77_BLAS_MANGLE(dpttrf, DPTTRF)

namespace ROL { 
  extern "C" {
    void PREFIX DGEEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info);
    void PREFIX DLATRS_F77(Teuchos_fcd UPLO, Teuchos_fcd TRANS, Teuchos_fcd DIAG, Teuchos_fcd NORMIN, const int* N, double* A, const int* LDA, double* X, double* SCALE, double* CNORM, int* INFO);
    void PREFIX DGTTRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* dl, const double* d, const double* du, const double* du2, const int* ipiv, double* x , const int* ldx, int* info);
    void PREFIX DGTTRF_F77(const int* n, double* dl, double* d, double* du, double* du2, int* ipiv, int* info);
    void PREFIX DSTEQR_F77(Teuchos_fcd, const int* n, double* D, double* E, double* Z, const int* ldz, double* work, int* info);
    void PREFIX DTRTRS_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda, double* b, const int* ldb, int* info);
    void PREFIX DGETRF_F77(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
    void PREFIX DGETRI_F77(const int* n, double* a, const int* lda, const int* ipiv, double* work , const int* lwork, int* info);
    void PREFIX DGETRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda,const int* ipiv, double* x , const int* ldx, int* info);
    void PREFIX DGEQRF_F77(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info);
    void PREFIX DORGQR_F77(const int* m, const int* n, const int* k, double* a, const int* lda, const double* tau, double* work, const int* lwork, int* info);
    void PREFIX DPTTRS_F77(const int* n, const int* nrhs, const double* d, const double* e, double* x , const int* ldx, int* info);
    void PREFIX DPTTRF_F77(const int* n, double* d, double* e, int* info);
  }

  template<typename OrdinalType, typename Real>
  struct LAPACK { 
    /// \brief Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues and, optionally, the left and/or right eigenvectors.
    ///
    /// Real and imaginary parts of the eigenvalues are returned in
    /// separate arrays, WR for real and WI for complex.  The RWORK
    /// array is only referenced if ScalarType is complex.
    void GEEV(const char& JOBVL, const char& JOBVR, const OrdinalType& n, Real* A, const OrdinalType& lda, Real* WR, Real* WI, Real* VL, const OrdinalType& ldvl, Real* VR, const OrdinalType& ldvr, Real* WORK, const OrdinalType& lwork, Real* RWORK, OrdinalType* info) const;

    /// \brief Robustly solve a possibly singular triangular linear system.
    ///
    /// \note This routine is slower than the BLAS' TRSM, but can
    ///   detect possible singularity of A.
    void LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN,
        const OrdinalType& N, Real* A, const OrdinalType& LDA, Real* X, Real* SCALE, Real* CNORM, OrdinalType* INFO) const;

    //! Computes an LU factorization of a \c n by \c n tridiagonal matrix \c A using partial pivoting with row interchanges.
    void GTTRF(const OrdinalType& n, Real* dl, Real* d, Real* du, Real* du2, OrdinalType* IPIV, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B or \c A^H*X=B with a tridiagonal matrix \c A using the LU factorization computed by GTTRF.
    void GTTRS(const char& TRANS, const OrdinalType& n, const OrdinalType& nrhs, const Real* dl,
        const Real* d, const Real* du, const Real* du2, const OrdinalType* IPIV, Real* B,
        const OrdinalType& ldb, OrdinalType* info) const;

    //! Computes the eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal \c n by \c n matrix \c A using implicit QL/QR.  The eigenvectors can only be computed if \c A was reduced to tridiagonal form by SYTRD.
    void STEQR(const char& COMPZ, const OrdinalType& n, Real* D, Real* E, Real* Z, const OrdinalType& ldz, Real* WORK, OrdinalType* info) const;

    //! Solves a triangular linear system of the form \c A*X=B or \c A**T*X=B, where \c A is a triangular matrix.
    void TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const OrdinalType& n, const OrdinalType& nrhs, const Real* A, const OrdinalType& lda, Real* B, const OrdinalType& ldb, OrdinalType* info) const;

    //! Computes an LU factorization of a general \c m by \c n matrix \c A using partial pivoting with row interchanges.
    void GETRF(const OrdinalType& m, const OrdinalType& n, Real* A, const OrdinalType& lda, OrdinalType* IPIV, OrdinalType* info) const;

    //! Computes the inverse of a matrix \c A using the LU factorization computed by GETRF.
    void GETRI(const OrdinalType& n, Real* A, const OrdinalType& lda, const OrdinalType* IPIV, Real* WORK, const OrdinalType& lwork, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general \c n by \c n matrix \c A using the LU factorization computed by GETRF.
    void GETRS(const char& TRANS, const OrdinalType& n, const OrdinalType& nrhs, const Real* A, const OrdinalType& lda, const OrdinalType* IPIV, Real* B, const OrdinalType& ldb, OrdinalType* info) const;

    //! Computes a QR factorization of a general \c m by \c n matrix \c A.
    void GEQRF (const OrdinalType& m, const OrdinalType& n, Real* A, const OrdinalType& lda, Real* TAU, Real* WORK, const OrdinalType& lwork, OrdinalType* info) const;

    /// \brief Compute explicit Q factor from QR factorization (GEQRF) (real case).
    ///
    /// Generate the \c m by \c n matrix Q with orthonormal columns
    /// corresponding to the first \c n columns of a product of \c k
    /// elementary reflectors of order \c m, as returned by \c GEQRF.
    ///
    /// \note This method is not defined when Real is complex.
    /// Call \c UNGQR in that case.  ("OR" stands for "orthogonal";
    /// "UN" stands for "unitary.")
    void ORGQR(const OrdinalType& m, const OrdinalType& n, const OrdinalType& k, Real* A, const OrdinalType& lda, const Real* TAU, Real* WORK, const OrdinalType& lwork, OrdinalType* info) const;

    //! Computes the \c L*D*L' factorization of a Hermitian/symmetric positive definite tridiagonal matrix \c A.
    void PTTRF(const OrdinalType& n, Real* d, Real* e, OrdinalType* info) const;

    //! Solves a tridiagonal system \c A*X=B using the \L*D*L' factorization of \c A computed by PTTRF.
    void PTTRS(const OrdinalType& n, const OrdinalType& nrhs, const Real* d, const Real* e, Real* B, const OrdinalType& ldb, OrdinalType* info) const;

  };

  template<>
  struct LAPACK<int, double> { 

    void LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN,
                const int& N, double* A, const int& LDA, double* X, double* SCALE, double* CNORM,
                int* INFO) const
    {
      DLATRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), CHAR_MACRO(NORMIN),
          &N, A, &LDA, X, SCALE, CNORM, INFO);
    }

    void GTTRS(const char& TRANS, const int& n, const int& nrhs, const double* dl, const double* d,
               const double* du, const double* du2, const int* IPIV, double* B, const int& ldb, int* info) const
    {
      DGTTRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, dl, d, du, du2, IPIV, B, &ldb, info); 
    }

    void GTTRF(const int& n, double* dl, double* d, double* du, double* du2, int* IPIV, int* info) const
    {
      DGTTRF_F77(&n, dl, d, du, du2, IPIV, info);
    }

    void STEQR(const char& COMPZ, const int& n, double* D, double* E, double* Z, const int& ldz, double* WORK, int* info) const
    {
      DSTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info);
    }

    void TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const int& n, const int& nrhs, const double* A, const int& lda, double* B, const int& ldb, int* info) const
    {
      DTRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &n, &nrhs, A, &lda, B, &ldb, info);
    }
    void GETRF(const int& m, const int& n, double* A, const int& lda, int* IPIV, int* info) const
    {
      DGETRF_F77(&m, &n, A, &lda, IPIV, info); 
    }

    void GETRI(const int& n, double* A, const int& lda, const int* IPIV, double* WORK, const int& lwork, int* info) const
    {
      DGETRI_F77(&n, A, &lda, IPIV, WORK, &lwork, info); 
    }

    void GETRS(const char& TRANS, const int& n, const int& nrhs, const double* A, const int& lda, const int* IPIV, double* B, const int& ldb, int* info) const
    {
      DGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info); 
    }

    void GEQRF( const int& m, const int& n, double* A, const int& lda, double* TAU, double* WORK, const int& lwork, int* info) const
    {
      DGEQRF_F77(&m, &n, A, &lda, TAU, WORK, &lwork, info); 
    }

    void ORGQR(const int& m, const int& n, const int& k, double* A, const int& lda, const double* TAU, double* WORK, const int& lwork, int* info) const
    {
      DORGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info);
    }

    void PTTRS(const int& n, const int& nrhs, const double* d, const double* e, double* B, const int& ldb, int* info) const
    { 
      DPTTRS_F77(&n,&nrhs,d,e,B,&ldb,info);
    }
    void GEEV(const char& JOBVL, const char& JOBVR, const int& n, double* A, const int& lda, double* WR, double* WI, double* VL, const int& ldvl, double* VR, const int& ldvr, double* WORK, const int& lwork, int* info) const
    {
      DGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, info);
    }
    void PTTRF(const int& n, double* d, double* e, int* info) const
    {
      DPTTRF_F77(&n,d,e,info); 
    }
  };

}

#endif
