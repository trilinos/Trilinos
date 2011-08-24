//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Tsqr_Blas.hpp>
#include <complex>

// C doesn't allow 'extern "C"' declarations inside a class' member
// functions, so we have to list all the 'extern "C"' declarations up
// here.

extern "C" void F77_BLAS_MANGLE(dgemv, DGEMV) 
  (const char* const TRANS,
   const int* const M,
   const int* const N,
   const double* const ALPHA,
   const double A[],
   const int* const LDA,
   const double X[],
   const int* const INCX,
   const double* const BETA,
   double Y[],
   const int* const INCY);

extern "C" void F77_BLAS_MANGLE(sgemv, SGEMV)
  (const char* const TRANS,
   const int* const M,
   const int* const N,
   const float* const ALPHA,
   const float A[],
   const int* const LDA,
   const float X[],
   const int* const INCX,
   const float* const BETA,
   float Y[],
   const int* const INCY);

extern "C" void F77_BLAS_MANGLE(zgemv, ZGEMV)
  (const char* const TRANS,
   const int* const M,
   const int* const N,
   const std::complex<double>* const ALPHA,
   const std::complex<double> A[],
   const int* const LDA,
   const std::complex<double> X[],
   const int* const INCX,
   const std::complex<double>* const BETA,
   std::complex<double> Y[],
   const int* const INCY);

extern "C" void F77_BLAS_MANGLE(cgemv, CGEMV)
  (const char* const TRANS,
   const int* const M,
   const int* const N,
   const std::complex<float>* const ALPHA,
   const std::complex<float> A[],
   const int* const LDA,
   const std::complex<float> X[],
   const int* const INCX,
   const std::complex<float>* const BETA,
   std::complex<float> Y[],
   const int* const INCY);

extern "C" void F77_BLAS_MANGLE(dgemm, DGEMM)
  (const char* const TRANSA,
   const char* const TRANSB,
   const int* const M,
   const int* const N,
   const int* const K,
   const double* const ALPHA,
   const double A[],
   const int* const LDA,
   const double B[],
   const int* const LDB,
   const double* const BETA,
   double C[],
   const int* const LDC);

extern "C" void F77_BLAS_MANGLE(sgemm, SGEMM)
  (const char* const TRANSA,
   const char* const TRANSB,
   const int* const M,
   const int* const N,
   const int* const K,
   const float* const ALPHA,
   const float A[],
   const int* const LDA,
   const float B[],
   const int* const LDB,
   const float* const BETA,
   float C[],
   const int* const LDC);

extern "C" void F77_BLAS_MANGLE(zgemm, ZGEMM)
  (const char* const TRANSA,
   const char* const TRANSB,
   const int* const M,
   const int* const N,
   const int* const K,
   const std::complex<double>* const ALPHA,
   const std::complex<double> A[],
   const int* const LDA,
   const std::complex<double> B[],
   const int* const LDB,
   const std::complex<double>* const BETA,
   std::complex<double> C[],
   const int* const LDC);

extern "C" void F77_BLAS_MANGLE(cgemm, CGEMM)
  (const char* const TRANSA,
   const char* const TRANSB,
   const int* const M,
   const int* const N,
   const int* const K,
   const std::complex<float>* const ALPHA,
   const std::complex<float> A[],
   const int* const LDA,
   const std::complex<float> B[],
   const int* const LDB,
   const std::complex<float>* const BETA,
   std::complex<float> C[],
   const int* const LDC);

extern "C" void F77_BLAS_MANGLE(dger, DGER)
  (const int* const M,
   const int* const N,
   const double* const ALPHA,
   const double X[],
   const int* const INCX,
   const double Y[],
   const int* const INCY,
   double A[],
   const int* const LDA);

extern "C" void F77_BLAS_MANGLE(sger, SGER)
  (const int* const M,
   const int* const N,
   const float* const ALPHA,
   const float X[],
   const int* const INCX,
   const float Y[],
   const int* const INCY,
   float A[],
   const int* const LDA);

extern "C" void F77_BLAS_MANGLE(zgerc, ZGERC)
  (const int* const M,
   const int* const N,
   const std::complex<double>* const ALPHA,
   const std::complex<double> X[],
   const int* const INCX,
   const std::complex<double> Y[],
   const int* const INCY,
   std::complex<double> A[],
   const int* const LDA);

extern "C" void F77_BLAS_MANGLE(cgerc, CGERC)
  (const int* const M,
   const int* const N,
   const std::complex<float>* const ALPHA,
   const std::complex<float> X[],
   const int* const INCX,
   const std::complex<float> Y[],
   const int* const INCY,
   std::complex<float> A[],
   const int* const LDA);

extern "C" void F77_BLAS_MANGLE(dtrsm, DTRSM)
  (const char* const SIDE,
   const char* const UPLO,
   const char* const TRANSA,
   const char* const DIAG,
   const int* const M,
   const int* const N,
   const double* const ALPHA,
   const double A[],
   const int* const LDA,
   double B[],
   const int* const LDB);

extern "C" void F77_BLAS_MANGLE(strsm, STRSM)
  (const char* const SIDE,
   const char* const UPLO,
   const char* const TRANSA,
   const char* const DIAG,
   const int* const M,
   const int* const N,
   const float* const ALPHA,
   const float A[],
   const int* const LDA,
   float B[],
   const int* const LDB);

extern "C" void F77_BLAS_MANGLE(ztrsm, ZTRSM)
  (const char* const SIDE,
   const char* const UPLO,
   const char* const TRANSA,
   const char* const DIAG,
   const int* const M,
   const int* const N,
   const std::complex<double>* const ALPHA,
   const std::complex<double> A[],
   const int* const LDA,
   std::complex<double> B[],
   const int* const LDB);

extern "C" void F77_BLAS_MANGLE(ctrsm, CTRSM)
  (const char* const SIDE,
   const char* const UPLO,
   const char* const TRANSA,
   const char* const DIAG,
   const int* const M,
   const int* const N,
   const std::complex<float>* const ALPHA,
   const std::complex<float> A[],
   const int* const LDA,
   std::complex<float> B[],
   const int* const LDB);

namespace TSQR {

  template<>
  void
  BLAS<int, double>::
  GEMV (const char* const trans, 
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
    F77_BLAS_MANGLE(dgemv, DGEMV) 
      (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, float>::
  GEMV (const char* const trans, 
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
    F77_BLAS_MANGLE(sgemv, SGEMV) 
      (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::
  GEMV (const char* const trans, 
	const int m, 
	const int n,
	const std::complex<double> alpha,
	const std::complex<double> A[],
	const int lda,
	const std::complex<double> x[],
	const int incx,
	const std::complex<double> beta,
	std::complex<double> y[],
	const int incy)
  {
    F77_BLAS_MANGLE(zgemv, ZGEMV) 
      (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::
  GEMV (const char* const trans, 
	const int m, 
	const int n,
	const std::complex<float> alpha,
	const std::complex<float> A[],
	const int lda,
	const std::complex<float> x[],
	const int incx,
	const std::complex<float> beta,
	std::complex<float> y[],
	const int incy)
  {
    F77_BLAS_MANGLE(cgemv, CGEMV) 
      (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, double>::
  GEMM (const char* const transa,
	const char* const transb,
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
    F77_BLAS_MANGLE(dgemm, DGEMM) 
      (transa, transb, &m, &n, &k, &alpha,
       A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, float>::
  GEMM (const char* const transa,
	const char* const transb,
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
    F77_BLAS_MANGLE(sgemm, SGEMM) 
      (transa, transb, &m, &n, &k, &alpha,
       A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::
  GEMM (const char* const transa,
	const char* const transb,
	const int m,
	const int n,
	const int k,
	const std::complex<double> alpha,
	const std::complex<double> A[],
	const int lda,
	const std::complex<double> B[],
	const int ldb,
	const std::complex<double> beta,
	std::complex<double> C[],
	const int ldc)
  {
    F77_BLAS_MANGLE(zgemm, ZGEMM) 
      (transa, transb, &m, &n, &k, &alpha,
       A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::
  GEMM (const char* const transa,
	const char* const transb,
	const int m,
	const int n,
	const int k,
	const std::complex<float> alpha,
	const std::complex<float> A[],
	const int lda,
	const std::complex<float> B[],
	const int ldb,
	const std::complex<float> beta,
	std::complex<float> C[],
	const int ldc)
  {
    F77_BLAS_MANGLE(cgemm, CGEMM) 
      (transa, transb, &m, &n, &k, &alpha, 
       A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, double>::
  GER (const int m,
       const int n,
       const double alpha,
       const double x[],
       const int incx,
       const double y[],
       const int incy,
       double A[],
       const int lda)
  {
    F77_BLAS_MANGLE(dger, DGER) 
      (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, float>::
  GER (const int m,
       const int n,
       const float alpha,
       const float x[],
       const int incx,
       const float y[],
       const int incy,
       float A[],
       const int lda)
  {
    F77_BLAS_MANGLE(sger, SGER) 
      (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::
  GER (const int m,
       const int n,
       const std::complex<double> alpha,
       const std::complex<double> x[],
       const int incx,
       const std::complex<double> y[],
       const int incy,
       std::complex<double> A[],
       const int lda)
  {
    F77_BLAS_MANGLE(zgerc, ZGERC) 
      (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::
  GER (const int m,
       const int n,
       const std::complex<float> alpha,
       const std::complex<float> x[],
       const int incx,
       const std::complex<float> y[],
       const int incy,
       std::complex<float> A[],
       const int lda)
  {
    F77_BLAS_MANGLE(cgerc, CGERC) 
      (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, double >::
  TRSM (const char* const side,
	const char* const uplo,
	const char* const transa,
	const char* const diag,
	const int m,
	const int n,
	const double alpha,
	const double A[],
	const int lda,
	double B[],
	const int ldb)
  {
    F77_BLAS_MANGLE(dtrsm, DTRSM) 
      (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

  template<>
  void
  BLAS<int, float >::
  TRSM (const char* const side,
	const char* const uplo,
	const char* const transa,
	const char* const diag,
	const int m,
	const int n,
	const float alpha,
	const float A[],
	const int lda,
	float B[],
	const int ldb)
  {
    F77_BLAS_MANGLE(strsm, STRSM) 
      (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::
  TRSM (const char* const side,
	const char* const uplo,
	const char* const transa,
	const char* const diag,
	const int m,
	const int n,
	const std::complex<double> alpha,
	const std::complex<double> A[],
	const int lda,
	std::complex<double> B[],
	const int ldb)
  {
    F77_BLAS_MANGLE(ztrsm, ZTRSM) 
      (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::
  TRSM (const char* const side,
	const char* const uplo,
	const char* const transa,
	const char* const diag,
	const int m,
	const int n,
	const std::complex<float> alpha,
	const std::complex<float> A[],
	const int lda,
	std::complex<float> B[],
	const int ldb)
  {
    F77_BLAS_MANGLE(ctrsm, CTRSM) 
      (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

} // namespace TSQR
