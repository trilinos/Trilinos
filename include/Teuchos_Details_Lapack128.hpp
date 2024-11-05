// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DETAILS_LAPACK128_HPP
#define TEUCHOS_DETAILS_LAPACK128_HPP

/// \file Teuchos_Details_Lapack128.hpp
/// \brief Declaration and definition of Teuchos::Details::Lapack128,
///   a partial implementation of Teuchos::LAPACK for __float128.

#include "Teuchos_ConfigDefs.hpp"


#ifdef HAVE_TEUCHOSCORE_QUADMATH
namespace Teuchos {
namespace Details {

//! Partial implementation of Teuchos::LAPACK for __float128.
class Lapack128 {
public:
  /// \brief Compute the LU factorization with partial pivoting of
  ///   the matrix A.
  void
  GETRF (const int M, const int N, __float128 A[],
         const int LDA, int IPIV[], int* INFO) const;

  /// \brief Perform a series of row interchanges on the matrix A.
  ///
  /// Do one row interchange for each of rows K1 through K2 of A.
  ///
  /// \param N [in] Number of columns of the matrix A.
  /// \param A [in/out] 2-D column-major array of dimension (LDA,N).
  ///   On entry, the matrix of column dimension N to which the row
  ///   interchanges will be applied.  On exit, the permuted matrix.
  /// \param LDA [in] The leading dimension (stride) of the 2-D
  ///   column-major array A.
  /// \param K1 [in] Start row interchanges with IPIV[K1-1].
  /// \param K2 [in] Stop row interchanges with IPIV[K2-1].
  /// \param INCX [in] Increment between successive entries of IPIV.
  ///   If IPIV is negative, apply the pivots in reverse order.
  void
  LASWP (const int N, __float128 A[], const int LDA, const int K1,
         const int K2, const int IPIV[], const int INCX) const;

  /// \brief Solve the linear system(s) AX=B, using the result of
  ///   the LU factorization computed by GETRF (above).
  void
  GETRS (const char TRANS, const int N, const int NRHS,
         const __float128 A[], const int LDA, const int IPIV[],
         __float128 B[], const int LDB, int* INFO) const;

  /// \brief Compute the inverse in place of the matrix A, using the
  ///   results of GETRF.
  void
  GETRI (const int N, __float128 A[], const int LDA, int IPIV[],
         __float128 WORK[], const int LWORK, int* INFO) const;

  /// \brief Compute the hypotenuse \f$\sqrt{x^2 + y^2}\f$ in a way
  ///   that avoids unjustified overflow.
  __float128
  LAPY2 (const __float128& x, const __float128& y) const;

  //! Compute the Householder reflector of [alpha; x].
  void
  LARFG (const int N, __float128* const ALPHA,
         __float128 X[], const int INCX, __float128* const TAU) const;

  //! Apply the Householder reflector [tau; v] to the matrix C.
  void
  LARF (const char side,
        const int m,
        const int n,
        const __float128 v[],
        const int incv,
        const __float128 tau,
        __float128 C[],
        const int ldc,
        __float128 work[]) const;

  //! BLAS 2 version of ORMQR; known workspace size.
  void
  ORM2R (const char side, const char trans,
         const int m, const int n, const int k,
         const __float128 A[], const int lda,
         const __float128* const tau,
         __float128 C[], const int ldc,
         __float128 work[], int* const info) const;

  //! BLAS 2 QR factorization of A.
  void
  GEQR2 (const int M,
         const int N,
         __float128 A[],
         const int LDA,
         __float128 TAU[],
         __float128 WORK[],
         int* const INFO) const;

  //! QR factorization of A.
  void
  GEQRF (const int M,
         const int N,
         __float128 A[],
         const int LDA,
         __float128 TAU[],
         __float128 WORK[],
         const int LWORK,
         int* const INFO) const;

  //! Assemble explicit Q factor from results of GEQRF (above).
  void
  ORGQR (const int M,
         const int N,
         const int K,
         __float128 A[],
         const int LDA,
         const __float128 TAU[],
         __float128 WORK[],
         const int LWORK,
         int* const INFO) const;

  //! Assemble explicit Q factor from results of GEQRF (above).
  void
  UNGQR (const int M,
         const int N,
         const int K,
         __float128 A[],
         const int LDA,
         const __float128 TAU[],
         __float128 WORK[],
         const int LWORK,
         int* const INFO) const;

  //! Scale the matrix A by the real scalar cto/cfrom.
  void
  LASCL (const char TYPE,
         const int kl,
         const int ku,
         const __float128 cfrom,
         const __float128 cto,
         const int m,
         const int n,
         __float128* A,
         const int lda,
         int* info) const;

  //! Compute LU factorization of the banded matrix A.
  void
  GBTRF (const int m,
         const int n,
         const int kl,
         const int ku,
         __float128* A,
         const int lda,
         int* IPIV,
         int* info) const;

  //! Solve linear system(s) using results of GBTRF (above).
  void
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
         int* info) const;
};

} // namespace Details
} // namespace Teuchos
#endif // HAVE_TEUCHOSCORE_QUADMATH

#endif // TEUCHOS_DETAILS_LAPACK128_HPP
