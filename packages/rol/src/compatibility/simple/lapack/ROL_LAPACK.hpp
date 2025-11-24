// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LAPACK_H
#define ROL_LAPACK_H

/** \class ROL::LAPACK
  \brief Provides interface to Lapack
  */

/* A) Define PREFIX and ROL_fcd based on platform. */


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
#  define ROL_fcd const char *, unsigned int
#elif defined(INTEL_MKL)
#  define PREFIX
#  define ROL_fcd const char *
#else
#  define PREFIX
#  define ROL_fcd const char *
#endif

// in progress - added for EQUED which is a modified char *, not const
#if defined(INTEL_CXML)
#  define PREFIX __stdcall
#  define ROL_nonconst_fcd char *, unsigned int // Need to evaluate unsigned int - CXML deprecated
#elif defined(INTEL_MKL)
#  define PREFIX
#  define ROL_nonconst_fcd char *
#else
#  define PREFIX
#  define ROL_nonconst_fcd char *
#endif

#define DGEEV_F77   F77_BLAS_MANGLE(dgeev,  DGEEV )
#define DGELS_F77   F77_BLAS_MANGLE(dgels,  DGELS )
#define DGELSS_F77  F77_BLAS_MANGLE(dgelss, DGELSS)
#define DGESVD_F77  F77_BLAS_MANGLE(dgesvd, DGESVD)
#define DGGEV_F77   F77_BLAS_MANGLE(dggev,  DGGEV)
#define DPOTRF_F77  F77_BLAS_MANGLE(dpotrf, DPOTRF)
#define DPOTRS_F77  F77_BLAS_MANGLE(dpotrs, DPOTRS)
#define DPORFS_F77  F77_BLAS_MANGLE(dporfs, DPORFS)
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
    void PREFIX DGEEV_F77(ROL_fcd, ROL_fcd, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info);
    void PREFIX DGELS_F77(ROL_fcd ch, const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb, double* work, const int* lwork, int* info);
    void PREFIX DGELSS_F77(const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb, double* s, const double* rcond, int* rank, double* work, const int* lwork, int* info);
    void PREFIX DGESVD_F77(ROL_fcd, ROL_fcd, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* v, const int* ldv, double* work, const int* lwork, int* info);
    void PREFIX DGGEV_F77(ROL_fcd, ROL_fcd, const int *n, double *A, const int *lda, double *B, const int *ldb, double *alphar, double *alphai, double *beta, double *vl, const int *ldvl, double *vr, const int *ldvr, double *work, const int *lwork, int *info);
    void PREFIX DPOTRF_F77(ROL_fcd, const int* n, double* a, const int* lda, int* info);
    void PREFIX DPOTRS_F77(ROL_fcd, const int* n, const int* nrhs, const double* a, const int* lda, double*x , const int* ldx, int* info);
    void PREFIX DPORFS_F77(ROL_fcd, const int* n, const int* nrhs, const double* a, const int* lda, const double* af, const int* ldaf, const double* b, const int* ldb, double* x, const int* ldx, double* ferr, double* berr, double* work, int* iwork, int* info);
    void PREFIX DLATRS_F77(ROL_fcd UPLO, ROL_fcd TRANS, ROL_fcd DIAG, ROL_fcd NORMIN, const int* N, double* A, const int* LDA, double* X, double* SCALE, double* CNORM, int* INFO);
    void PREFIX DGTTRS_F77(ROL_fcd, const int* n, const int* nrhs, const double* dl, const double* d, const double* du, const double* du2, const int* ipiv, double* x , const int* ldx, int* info);
    void PREFIX DGTTRF_F77(const int* n, double* dl, double* d, double* du, double* du2, int* ipiv, int* info);
    void PREFIX DSTEQR_F77(ROL_fcd, const int* n, double* D, double* E, double* Z, const int* ldz, double* work, int* info);
    void PREFIX DTRTRS_F77(ROL_fcd, ROL_fcd, ROL_fcd, const int* n, const int* nrhs, const double* a, const int* lda, double* b, const int* ldb, int* info);
    void PREFIX DGETRF_F77(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
    void PREFIX DGETRI_F77(const int* n, double* a, const int* lda, const int* ipiv, double* work , const int* lwork, int* info);
    void PREFIX DGETRS_F77(ROL_fcd, const int* n, const int* nrhs, const double* a, const int* lda,const int* ipiv, double* x , const int* ldx, int* info);
    void PREFIX DGEQRF_F77(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info);
    void PREFIX DORGQR_F77(const int* m, const int* n, const int* k, double* a, const int* lda, const double* tau, double* work, const int* lwork, int* info);
    void PREFIX DPTTRS_F77(const int* n, const int* nrhs, const double* d, const double* e, double* x , const int* ldx, int* info);
    void PREFIX DPTTRF_F77(const int* n, double* d, double* e, int* info);
  }

  template<typename Index, typename Real>
  struct LAPACK {
    /// \brief Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues and, optionally, the left and/or right eigenvectors.
    ///
    /// Real and imaginary parts of the eigenvalues are returned in
    /// separate arrays, WR for real and WI for complex.  The RWORK
    /// array is only referenced if ScalarType is complex.
    void GEEV(const char& JOBVL, const char& JOBVR, const Index& n, Real* A, const Index& lda, Real* WR, Real* WI, Real* VL, const Index& ldvl, Real* VR, const Index& ldvr, Real* WORK, const Index& lwork, Real* RWORK, Index* info) const;

    //\brief Solves an over/underdetermined real \c m by \c n linear system \c A using QR or LQ factorization of A.
    void GELS(const char& TRANS, const Index& m, const Index& n, const Index& nrhs, Real* A, const Index& lda, Real* B, const Index& ldb, Real* WORK, const Index& lwork, Index* info) const;

    void GELSS(const Index& m, const Index& n, const Index& nrhs, Real* A, const Index& lda, Real* B, const Index& ldb, Real* S, const Real& rcond, Index* rank, Real* WORK, const Index& lwork, Index* info) const;

    void GESVD(const char& JOBU, const char& JOBVT, const Index& m, const Index& n, Real* A, const Index& lda, Real* S, Real* U, const Index& ldu, Real* V, const Index& ldv, Real* WORK, const Index& lwork, Real* /* RWORK */, Index* info) const;

    /*! Computes for a pair of \c n by \c n nonsymmetric matrices (\c A,\c B) the generalized eigenvalues, and optionally, the left and/or right generalized eigenvectors.
       \note (This is the function is only defined for \c Real = \c float or \c double.)
    */
    void GGEV(const char& JOBVL, const char& JOBVR, const Index& n, Real* A, const Index& lda, Real* B, const Index& ldb, Real *ALPHAR, Real *ALPHAI, Real* BETA, Real* VL, const Index& ldvl, Real* VR, const Index& ldvr, Real* WORK, const Index& lwork, Index* info) const;

    //! Computes Cholesky factorization of a real symmetric positive definite matrix \c A.
    void POTRF(const char& UPLO, const Index& n, Real* A, const Index& lda, Index* info) const;

    //! Solves a system of linear equations \c A*X=B, where \c A is a symmetric positive definite matrix factored by POTRF and the \c nrhs solutions are returned in \c B.
    void POTRS(const char& UPLO, const Index& n, const Index& nrhs, const Real* A, const Index& lda, Real* B, const Index& ldb, Index* info) const;

    //! Improves the computed solution to a system of linear equations when the coefficient matrix is symmetric positive definite, and provides error bounds and backward error estimates for the solution.
    void PORFS(const char& UPLO, const Index& n, const Index& nrhs, const Real* A, const Index& lda, const Real* AF, const Index& ldaf, const Real* B, const Index& ldb, Real* X, const Index& ldx, Real* FERR, Real* BERR, Real* WORK, Index* IWORK, Index* info) const;

    /// \brief Robustly solve a possibly singular triangular linear system.
    ///
    /// \note This routine is slower than the BLAS' TRSM, but can
    ///   detect possible singularity of A.
    void LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN,
        const Index& N, Real* A, const Index& LDA, Real* X, Real* SCALE, Real* CNORM, Index* INFO) const;

    //! Computes an LU factorization of a \c n by \c n tridiagonal matrix \c A using partial pivoting with row interchanges.
    void GTTRF(const Index& n, Real* dl, Real* d, Real* du, Real* du2, Index* IPIV, Index* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B or \c A^H*X=B with a tridiagonal matrix \c A using the LU factorization computed by GTTRF.
    void GTTRS(const char& TRANS, const Index& n, const Index& nrhs, const Real* dl,
        const Real* d, const Real* du, const Real* du2, const Index* IPIV, Real* B,
        const Index& ldb, Index* info) const;

    //! Computes the eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal \c n by \c n matrix \c A using implicit QL/QR.  The eigenvectors can only be computed if \c A was reduced to tridiagonal form by SYTRD.
    void STEQR(const char& COMPZ, const Index& n, Real* D, Real* E, Real* Z, const Index& ldz, Real* WORK, Index* info) const;

    //! Solves a triangular linear system of the form \c A*X=B or \c A**T*X=B, where \c A is a triangular matrix.
    void TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const Index& n, const Index& nrhs, const Real* A, const Index& lda, Real* B, const Index& ldb, Index* info) const;

    //! Computes an LU factorization of a general \c m by \c n matrix \c A using partial pivoting with row interchanges.
    void GETRF(const Index& m, const Index& n, Real* A, const Index& lda, Index* IPIV, Index* info) const;

    //! Computes the inverse of a matrix \c A using the LU factorization computed by GETRF.
    void GETRI(const Index& n, Real* A, const Index& lda, const Index* IPIV, Real* WORK, const Index& lwork, Index* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general \c n by \c n matrix \c A using the LU factorization computed by GETRF.
    void GETRS(const char& TRANS, const Index& n, const Index& nrhs, const Real* A, const Index& lda, const Index* IPIV, Real* B, const Index& ldb, Index* info) const;

    //! Computes a QR factorization of a general \c m by \c n matrix \c A.
    void GEQRF (const Index& m, const Index& n, Real* A, const Index& lda, Real* TAU, Real* WORK, const Index& lwork, Index* info) const;

    /// \brief Compute explicit Q factor from QR factorization (GEQRF) (real case).
    ///
    /// Generate the \c m by \c n matrix Q with orthonormal columns
    /// corresponding to the first \c n columns of a product of \c k
    /// elementary reflectors of order \c m, as returned by \c GEQRF.
    ///
    /// \note This method is not defined when Real is complex.
    /// Call \c UNGQR in that case.  ("OR" stands for "orthogonal";
    /// "UN" stands for "unitary.")
    void ORGQR(const Index& m, const Index& n, const Index& k, Real* A, const Index& lda, const Real* TAU, Real* WORK, const Index& lwork, Index* info) const;

    //! Computes the \c L*D*L' factorization of a Hermitian/symmetric positive definite tridiagonal matrix \c A.
    void PTTRF(const Index& n, Real* d, Real* e, Index* info) const;

    //! Solves a tridiagonal system \c A*X=B using the \L*D*L' factorization of \c A computed by PTTRF.
    void PTTRS(const Index& n, const Index& nrhs, const Real* d, const Real* e, Real* B, const Index& ldb, Index* info) const;
  };

  template<>
  struct LAPACK<int, double> {

    void LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN,
                const int& N, double* A, const int& LDA, double* X, double* SCALE, double* CNORM,
                int* INFO) const {
      DLATRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), CHAR_MACRO(NORMIN),
          &N, A, &LDA, X, SCALE, CNORM, INFO);
    }

    void GTTRS(const char& TRANS, const int& n, const int& nrhs, const double* dl, const double* d,
               const double* du, const double* du2, const int* IPIV, double* B, const int& ldb, int* info) const {
      DGTTRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, dl, d, du, du2, IPIV, B, &ldb, info);
    }

    void GTTRF(const int& n, double* dl, double* d, double* du, double* du2, int* IPIV, int* info) const {
      DGTTRF_F77(&n, dl, d, du, du2, IPIV, info);
    }

    void STEQR(const char& COMPZ, const int& n, double* D, double* E, double* Z, const int& ldz, double* WORK, int* info) const {
      DSTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info);
    }

    void TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const int& n, const int& nrhs, const double* A, const int& lda, double* B, const int& ldb, int* info) const {
      DTRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &n, &nrhs, A, &lda, B, &ldb, info);
    }

    void GETRF(const int& m, const int& n, double* A, const int& lda, int* IPIV, int* info) const {
      DGETRF_F77(&m, &n, A, &lda, IPIV, info);
    }

    void GETRI(const int& n, double* A, const int& lda, const int* IPIV, double* WORK, const int& lwork, int* info) const {
      DGETRI_F77(&n, A, &lda, IPIV, WORK, &lwork, info);
    }

    void GETRS(const char& TRANS, const int& n, const int& nrhs, const double* A, const int& lda, const int* IPIV, double* B, const int& ldb, int* info) const {
      DGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info);
    }

    void GEQRF( const int& m, const int& n, double* A, const int& lda, double* TAU, double* WORK, const int& lwork, int* info) const {
      DGEQRF_F77(&m, &n, A, &lda, TAU, WORK, &lwork, info);
    }

    void ORGQR(const int& m, const int& n, const int& k, double* A, const int& lda, const double* TAU, double* WORK, const int& lwork, int* info) const {
      DORGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info);
    }

    void PTTRS(const int& n, const int& nrhs, const double* d, const double* e, double* B, const int& ldb, int* info) const {
      DPTTRS_F77(&n,&nrhs,d,e,B,&ldb,info);
    }

    void GEEV(const char& JOBVL, const char& JOBVR, const int& n, double* A, const int& lda, double* WR, double* WI, double* VL, const int& ldvl, double* VR, const int& ldvr, double* WORK, const int& lwork, int* info) const {
      DGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, info);
    }

    void GELS(const char& TRANS, const int& m, const int& n, const int& nrhs, double* A, const int& lda, double* B, const int& ldb, double* WORK, const int& lwork, int* info) const {
      DGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, A, &lda, B, &ldb, WORK, &lwork, info);
    }

    void GELSS(const int& m, const int& n, const int& nrhs, double* A, const int& lda, double* B, const int& ldb, double* S, const double& rcond, int* rank, double* WORK, const int& lwork, int* info) const {
      DGELSS_F77(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank, WORK, &lwork, info);
    }

    void GESVD(const char& JOBU, const char& JOBVT, const int& m, const int& n, double* A, const int& lda, double* S, double* U, const int& ldu, double* V, const int& ldv, double* WORK, const int& lwork, double* /* RWORK */, int* info) const {
      DGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &m, &n, A, &lda, S, U, &ldu, V, &ldv, WORK, &lwork, info);
    }

    void GGEV(const char& JOBVL, const char& JOBVR, const int& n, double* A, const int& lda, double* B, const int& ldb, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, const int& ldvl, double* VR, const int& ldvr, double* WORK, const int& lwork, int* info) const {
      DGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, B, &ldb, ALPHAR, ALPHAI, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, info);
    }

    void POTRF(const char& UPLO, const int& n, double* A, const int& lda, int* info) const {
      DPOTRF_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
    }

    void POTRS(const char& UPLO, const int& n, const int& nrhs, const double* A, const int& lda, double* B, const int& ldb, int* info) const {
      DPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
    }

    void PORFS(const char& UPLO, const int& n, const int& nrhs, const double* A, const int& lda, const double* AF, const int& ldaf, const double* B, const int& ldb, double* X, const int& ldx, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const {
      DPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info);
    }

    void PTTRF(const int& n, double* d, double* e, int* info) const {
      DPTTRF_F77(&n,d,e,info);
    }
  };

}

#endif
