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

#include "Tpetra_Details_libGemv.hpp"
#include "KokkosKernels_config.h"
#include <stdexcept>

// TpetraCore_config.h (included in Tpetra_Details_libGemv.hpp)
// defines the TPETRACORE_F77_BLAS_MANGLE macro.  First argument is
// the lower-case version of the BLAS function name, and second
// argument is the upper-case version of the same name.  The macro
// handles mangling for calling the Fortran 77 function from C.  We
// then make an extern "C" declaration here for each BLAS function
// that we want to use.

#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT

#define TPETRACORE_CGEMV TPETRACORE_F77_BLAS_MANGLE(cgemv,CGEMV)

extern "C" void
TPETRACORE_CGEMV
 (const char* trans,
  const int* m,
  const int* n,
  const void* alpha,
  const void* a,
  const int* lda,
  const void* x,
  const int* incx,
  const void* beta,
  void* y,
  const int* incy);

#endif // HAVE_TPETRA_INST_COMPLEX_FLOAT

#define TPETRACORE_DGEMV TPETRACORE_F77_BLAS_MANGLE(dgemv,DGEMV)

extern "C" void
TPETRACORE_DGEMV
 (const char* trans,
  const int* m,
  const int* n,
  const double* alpha,
  const double* a,
  const int* lda,
  const double* x,
  const int* incx,
  const double* beta,
  double* y,
  const int* incy);

#define TPETRACORE_SGEMV TPETRACORE_F77_BLAS_MANGLE(sgemv,SGEMV)

extern "C" void
TPETRACORE_SGEMV
  (const char* trans,
   const int* m,
   const int* n,
   const float* alpha,
   const float* a,
   const int* lda,
   const float* x,
   const int* incx,
   const float* beta,
   float* y,
   const int* incy);

#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE

#define TPETRACORE_ZGEMV TPETRACORE_F77_BLAS_MANGLE(zgemv,ZGEMV)

extern "C" void
TPETRACORE_ZGEMV
  (const char* trans,
   const int* m,
   const int* n,
   const void* alpha,
   const void* a,
   const int* lda,
   const void* x,
   const int* incx,
   const void* beta,
   void* y,
   const int* incy);

#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Lib {
namespace Impl {

void
cgemv (const char trans,
       const int m,
       const int n,
       const ::Kokkos::complex<float>& alpha,
       const ::Kokkos::complex<float> A[],
       const int lda,
       const ::Kokkos::complex<float> x[],
       const int incx,
       const ::Kokkos::complex<float>& beta,
       ::Kokkos::complex<float> y[],
       const int incy)
{
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  TPETRACORE_CGEMV (&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
#else
  throw std::runtime_error ("Tpetra::Details::Blas::Lib::Impl::cgemv: "
                            "You must configure Tpetra with complex<float> support.");
#endif // HAVE_TPETRA_INST_COMPLEX_FLOAT
}

void
dgemv (const char trans,
       const int m,
       const int n,
       const double alpha,
       const double A[],
       const int lda,
       const double x[],
       const int incx,
       const double beta,
       double y[],
       const int incy)
{
  TPETRACORE_DGEMV (&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

void
sgemv (const char trans,
       const int m,
       const int n,
       const float alpha,
       const float A[],
       const int lda,
       const float x[],
       const int incx,
       const float beta,
       float y[],
       const int incy)
{
  TPETRACORE_SGEMV (&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

void
zgemv (const char trans,
       const int m,
       const int n,
       const ::Kokkos::complex<double>& alpha,
       const ::Kokkos::complex<double> A[],
       const int lda,
       const ::Kokkos::complex<double> x[],
       const int incx,
       const ::Kokkos::complex<double>& beta,
       ::Kokkos::complex<double> y[],
       const int incy)
{
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  TPETRACORE_ZGEMV (&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
#else
  throw std::runtime_error ("Tpetra::Details::Blas::Lib::Impl::zgemv: "
                            "You must configure Tpetra with complex<double> support.");
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE
}

} // namespace Impl
} // namespace Lib
} // namespace Blas
} // namespace Details
} // namespace Tpetra
