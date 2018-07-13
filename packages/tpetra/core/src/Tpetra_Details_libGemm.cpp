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

#include "Tpetra_Details_libGemm.hpp"
#include "KokkosKernels_config.h"
#include <stdexcept>

// TpetraCore_config.h (included in Tpetra_Details_libGemm.hpp)
// defines the TPETRACORE_F77_BLAS_MANGLE macro.  First argument is
// the lower-case version of the BLAS function name, and second
// argument is the upper-case version of the same name.  The macro
// handles mangling for calling the Fortran 77 function from C.  We
// then make an extern "C" declaration here for each BLAS function
// that we want to use.

#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT

#define TPETRACORE_CGEMM TPETRACORE_F77_BLAS_MANGLE(cgemm,CGEMM)

extern "C" void
TPETRACORE_CGEMM
 (const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const void* alpha, const void* a, const int* lda,
  const void* b, const int* ldb,
  const void* beta, void* c, const int* ldc);

#endif // HAVE_TPETRA_INST_COMPLEX_FLOAT

#define TPETRACORE_DGEMM TPETRACORE_F77_BLAS_MANGLE(dgemm,DGEMM)

extern "C" void
TPETRACORE_DGEMM
 (const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const double* alpha, const double* a, const int* lda,
  const double* b, const int* ldb,
  const double* beta, double* c, const int* ldc);

#define TPETRACORE_SGEMM TPETRACORE_F77_BLAS_MANGLE(sgemm,SGEMM)

extern "C" void
TPETRACORE_SGEMM
  (const char* transA, const char* transB,
   const int* m, const int* n, const int* k,
   const float* alpha, const float* a, const int* lda,
   const float* b, const int* ldb,
   const float* beta, float* c, const int* ldc);

#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE

#define TPETRACORE_ZGEMM TPETRACORE_F77_BLAS_MANGLE(zgemm,ZGEMM)

extern "C" void
TPETRACORE_ZGEMM
  (const char* transA, const char* transB,
   const int* m, const int* n, const int* k,
   const void* alpha, const void* a, const int* lda,
   const void* b, const int* ldb,
   const void* beta, void* c, const int* ldc);

#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Lib {
namespace Impl {

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
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  TPETRACORE_CGEMM (&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
#else
  throw std::runtime_error ("Tpetra::Details::Blas::Lib::Impl::cgemm: "
                            "You must configure Tpetra with complex<float> support.");
#endif // HAVE_TPETRA_INST_COMPLEX_FLOAT
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
  TPETRACORE_DGEMM (&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
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
  TPETRACORE_SGEMM (&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
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
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  TPETRACORE_ZGEMM (&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
#else
  throw std::runtime_error ("Tpetra::Details::Blas::Lib::Impl::zgemm: "
                            "You must configure Tpetra with complex<double> support.");
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE
}

} // namespace Impl
} // namespace Lib
} // namespace Blas
} // namespace Details
} // namespace Tpetra
