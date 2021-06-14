// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_DETAILS_LAPACKLONGDOUBLE_HPP
#define TEUCHOS_DETAILS_LAPACKLONGDOUBLE_HPP

/// \file Teuchos_Details_LapackLongDouble.hpp
/// \brief Declaration and definition of Teuchos::Details::LapackLongDouble,
///   a partial implementation of Teuchos::LAPACK for long double.

#include "Teuchos_ConfigDefs.hpp"


#ifdef HAVE_TEUCHOS_LONG_DOUBLE
namespace Teuchos {
namespace Details {

//! Partial implementation of Teuchos::LAPACK for long double.
class LapackLongDouble {
public:
  /// \brief Compute the LU factorization with partial pivoting of
  ///   the matrix A.
  void
  GETRF (const int M, const int N, long double A[],
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
  LASWP (const int N, long double A[], const int LDA, const int K1,
         const int K2, const int IPIV[], const int INCX) const;

  /// \brief Solve the linear system(s) AX=B, using the result of
  ///   the LU factorization computed by GETRF (above).
  void
  GETRS (const char TRANS, const int N, const int NRHS,
         const long double A[], const int LDA, const int IPIV[],
         long double B[], const int LDB, int* INFO) const;

  /// \brief Compute the inverse in place of the matrix A, using the
  ///   results of GETRF.
  void
  GETRI (const int N, long double A[], const int LDA, int IPIV[],
         long double WORK[], const int LWORK, int* INFO) const;

  /// \brief Compute the hypotenuse \f$\sqrt{x^2 + y^2}\f$ in a way
  ///   that avoids unjustified overflow.
  long double
  LAPY2 (const long double& x, const long double& y) const;

  //! Compute the Householder reflector of [alpha; x].
  void
  LARFG (const int N, long double* const ALPHA,
         long double X[], const int INCX, long double* const TAU) const;

  //! Apply the Householder reflector [tau; v] to the matrix C.
  void
  LARF (const char side,
        const int m,
        const int n,
        const long double v[],
        const int incv,
        const long double tau,
        long double C[],
        const int ldc,
        long double work[]) const;

  //! BLAS 2 version of ORMQR; known workspace size.
  void
  ORM2R (const char side, const char trans,
         const int m, const int n, const int k,
         const long double A[], const int lda,
         const long double* const tau,
         long double C[], const int ldc,
         long double work[], int* const info) const;

  //! BLAS 2 QR factorization of A.
  void
  GEQR2 (const int M,
         const int N,
         long double A[],
         const int LDA,
         long double TAU[],
         long double WORK[],
         int* const INFO) const;

  //! QR factorization of A.
  void
  GEQRF (const int M,
         const int N,
         long double A[],
         const int LDA,
         long double TAU[],
         long double WORK[],
         const int LWORK,
         int* const INFO) const;

  //! Assemble explicit Q factor from results of GEQRF (above).
  void
  ORGQR (const int M,
         const int N,
         const int K,
         long double A[],
         const int LDA,
         const long double TAU[],
         long double WORK[],
         const int LWORK,
         int* const INFO) const;

  //! Assemble explicit Q factor from results of GEQRF (above).
  void
  UNGQR (const int M,
         const int N,
         const int K,
         long double A[],
         const int LDA,
         const long double TAU[],
         long double WORK[],
         const int LWORK,
         int* const INFO) const;

  //! Scale the matrix A by the real scalar cto/cfrom.
  void
  LASCL (const char TYPE,
         const int kl,
         const int ku,
         const long double cfrom,
         const long double cto,
         const int m,
         const int n,
         long double* A,
         const int lda,
         int* info) const;

  //! Compute LU factorization of the banded matrix A.
  void
  GBTRF (const int m,
         const int n,
         const int kl,
         const int ku,
         long double* A,
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
         const long double* A,
         const int lda,
         const int* IPIV,
         long double* B,
         const int ldb,
         int* info) const;
};

} // namespace Details
} // namespace Teuchos
#endif // HAVE_TEUCHOS_LONG_DOUBLE

#endif // TEUCHOS_DETAILS_LAPACKLONGDOUBLE_HPP
