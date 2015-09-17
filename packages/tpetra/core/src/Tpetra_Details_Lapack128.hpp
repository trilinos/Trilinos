// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_LAPACK128_HPP
#define TPETRA_DETAILS_LAPACK128_HPP

/// \file Tpetra_Details_Lapack128.hpp
/// \brief Declaration and definition of Tpetra::Details::Lapack128,
///   a partial implementation of Teuchos::LAPACK for __float128.

#include "Tpetra_ConfigDefs.hpp"


#ifdef HAVE_TPETRA_INST_FLOAT128
namespace Tpetra {
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
};

} // namespace Details
} // namespace Tpetra
#endif // HAVE_TPETRA_INST_FLOAT128

#endif // TPETRA_DETAILS_LAPACK128_HPP
