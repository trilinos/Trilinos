/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#include "Tpetra_Details_mklGemm.hpp"
#include "KokkosKernels_config.h"
#ifdef HAVE_KOKKOSKERNELS_MKL
#  include <mkl.h>
#endif // HAVE_KOKKOSKERNELS_MKL
#include <sstream>
#include <stdexcept>

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Mkl {
namespace Impl {

#ifdef HAVE_KOKKOSKERNELS_MKL
namespace { // (anonymous)

CBLAS_TRANSPOSE
transCharToMklEnum (const char trans)
{
  if (trans == 'T' || trans == 't') {
    return CblasTrans;
  }
  else if (trans == 'C' || trans == 'c' ||
           trans == 'H' || trans == 'h') {
    return CblasConjTrans;
  }
  else {
    return CblasNoTrans;
  }
}

} // namespace (anonymous)
#endif // HAVE_KOKKOSKERNELS_MKL

void
cgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const ::Kokkos::complex<float>& alpha,
       const ::Kokkos::complex<float> A[],
       const int lda,
       const ::Kokkos::complex<float> B[],
       const int ldb,
       const ::Kokkos::complex<float>& beta,
       ::Kokkos::complex<float> C[],
       const int ldc)
{
#ifdef HAVE_KOKKOSKERNELS_MKL
  const CBLAS_TRANSPOSE etransA = transCharToMklEnum (transA);
  const CBLAS_TRANSPOSE etransB = transCharToMklEnum (transB);
  // MKL's complex interface, unlike real interface, takes alpha and
  // beta by pointer.
  ::cblas_cgemm (CblasColMajor, etransA, etransB,
                 m, n, k,
                 &alpha, A, lda,
                 B, ldb,
                 &beta, C, ldc);
#else // NOT HAVE_KOKKOSKERNELS_MKL
  throw std::runtime_error ("You must enable MKL in your Trilinos build in "
                            "order to invoke MKL functions in Tpetra.");
#endif // NOT HAVE_KOKKOSKERNELS_MKL
}

void
dgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const double alpha,
       const double A[],
       const int lda,
       const double B[],
       const int ldb,
       const double beta,
       double C[],
       const int ldc)
{
#ifdef HAVE_KOKKOSKERNELS_MKL
  const CBLAS_TRANSPOSE etransA = transCharToMklEnum (transA);
  const CBLAS_TRANSPOSE etransB = transCharToMklEnum (transB);
  ::cblas_dgemm (CblasColMajor, etransA, etransB,
                 m, n, k,
                 alpha, A, lda,
                 B, ldb,
                 beta, C, ldc);
#else // NOT HAVE_KOKKOSKERNELS_MKL
  throw std::runtime_error ("You must enable MKL in your Trilinos build in "
                            "order to invoke MKL functions in Tpetra.");
#endif // NOT HAVE_KOKKOSKERNELS_MKL
}

void
sgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const float alpha,
       const float A[],
       const int lda,
       const float B[],
       const int ldb,
       const float beta,
       float C[],
       const int ldc)
{
#ifdef HAVE_KOKKOSKERNELS_MKL
  const CBLAS_TRANSPOSE etransA = transCharToMklEnum (transA);
  const CBLAS_TRANSPOSE etransB = transCharToMklEnum (transB);
  ::cblas_sgemm (CblasColMajor, etransA, etransB,
                 m, n, k,
                 alpha, A, lda,
                 B, ldb,
                 beta, C, ldc);
#else // NOT HAVE_KOKKOSKERNELS_MKL
  throw std::runtime_error ("You must enable MKL in your Trilinos build in "
                            "order to invoke MKL functions in Tpetra.");
#endif // NOT HAVE_KOKKOSKERNELS_MKL
}

void
zgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const ::Kokkos::complex<double>& alpha,
       const ::Kokkos::complex<double> A[],
       const int lda,
       const ::Kokkos::complex<double> B[],
       const int ldb,
       const ::Kokkos::complex<double>& beta,
       ::Kokkos::complex<double> C[],
       const int ldc)
{
#ifdef HAVE_KOKKOSKERNELS_MKL
  const CBLAS_TRANSPOSE etransA = transCharToMklEnum (transA);
  const CBLAS_TRANSPOSE etransB = transCharToMklEnum (transB);
  // MKL's complex interface, unlike real interface, takes alpha and
  // beta by pointer.
  ::cblas_zgemm (CblasColMajor, etransA, etransB,
                 m, n, k,
                 &alpha, A, lda,
                 B, ldb,
                 &beta, C, ldc);
#else // NOT HAVE_KOKKOSKERNELS_MKL
  throw std::runtime_error ("You must enable MKL in your Trilinos build in "
                            "order to invoke MKL functions in Tpetra.");
#endif // NOT HAVE_KOKKOSKERNELS_MKL
}

} // namespace Impl
} // namespace Mkl
} // namespace Blas
} // namespace Details
} // namespace Tpetra
