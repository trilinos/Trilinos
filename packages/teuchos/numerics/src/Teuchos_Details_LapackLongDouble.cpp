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

#include "Teuchos_Details_LapackLongDouble.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef HAVE_TEUCHOS_LONG_DOUBLE
namespace Teuchos {
namespace Details {

void
LapackLongDouble::
GETRF (const int M, const int N, long double A[],
       const int LDA, int IPIV[], int* INFO) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GETRF: Not implemented yet.");
}

void
LapackLongDouble::
LASWP (const int N, long double A[], const int LDA, const int K1,
       const int K2, const int IPIV[], const int INCX) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::LASWP: Not implemented yet.");
}

void
LapackLongDouble::
GETRI (const int /* N */, long double /* A */ [], const int /* LDA */,
       int /* IPIV */ [], long double /* WORK */ [], const int /* LWORK */,
       int* /* INFO */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GETRI: Not implemented yet.");
}


void
LapackLongDouble::
GETRS (const char TRANS, const int N, const int NRHS,
       const long double A[], const int LDA, const int IPIV[],
       long double B[], const int LDB, int* INFO) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GETRS: Not implemented yet.");
}

long double
LapackLongDouble::
LAPY2 (const long double& x, const long double& y) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::LAPY2: Not implemented yet.");
}

void
LapackLongDouble::
ORM2R (const char side, const char trans,
       const int m, const int n, const int k,
       const long double A[], const int lda,
       const long double* const tau,
       long double C[], const int ldc,
       long double work[], int* const info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::ORM2R: Not implemented yet.");
}

namespace { // (anonymous)

  int
  ILADLC (const int m, const int n, const long double A[], const int lda)
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Teuchos::LAPACK<int, long double>::ILADLC: Not implemented yet.");
    return 0;
  }

  int
  ILADLR (const int m, const int n, const long double A[], const int lda)
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Teuchos::LAPACK<int, long double>::ILADLR: Not implemented yet.");
    return 0; 
  }
} // namespace (anonymous)

void
LapackLongDouble::
LARF (const char side,
      const int m,
      const int n,
      const long double v[],
      const int incv,
      const long double tau,
      long double C[],
      const int ldc,
      long double work[]) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::LARF: Not implemented yet.");
}


void
LapackLongDouble::
LARFG (const int N, long double* const ALPHA,
       long double X[], const int INCX, long double* const TAU) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::LARFG: Not implemented yet.");
}

void
LapackLongDouble::
GEQR2 (const int /* M */,
       const int /* N */,
       long double /* A */ [],
       const int /* LDA */,
       long double /* TAU */ [],
       long double /* WORK */ [],
       int* const /* INFO */ ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GEQR2: Not implemented yet.");
}

void
LapackLongDouble::
GEQRF (const int M,
       const int N,
       long double A[],
       const int LDA,
       long double TAU[],
       long double WORK[],
       const int LWORK,
       int* const INFO) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GEQRF: Not implemented yet.");
}

void
LapackLongDouble::
ORGQR (const int /* M */,
       const int /* N */,
       const int /* K */,
       long double /* A */ [],
       const int /* LDA */,
       const long double /* TAU */ [],
       long double /* WORK */ [],
       const int /* LWORK */,
       int* const /* INFO */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GEQR2: Not implemented yet.");
}

void
LapackLongDouble::
UNGQR (const int /* M */,
       const int /* N */,
       const int /* K */,
       long double /* A */ [],
       const int /* LDA */,
       const long double /* TAU */ [],
       long double /* WORK */ [],
       const int /* LWORK */,
       int* const /* INFO */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GEQR2: Not implemented yet.");
}

void
LapackLongDouble::
LASCL (const char TYPE,
       const int kl,
       const int ku,
       const long double cfrom,
       const long double cto,
       const int m,
       const int n,
       long double* A,
       const int lda,
       int* info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::LASCL: Not implemented yet.");
}

void
LapackLongDouble::
GBTRF (const int m,
       const int n,
       const int kl,
       const int ku,
       long double* A,
       const int lda,
       int* IPIV,
       int* info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GBTRF: Not implemented yet.");
}

void
LapackLongDouble::
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
       int* info) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Teuchos::LAPACK<int, long double>::GBTRS: Not implemented yet.");
}

} // namespace Details
} // namespace Teuchos
#endif // HAVE_TEUCHOS_LONG_DOUBLE
