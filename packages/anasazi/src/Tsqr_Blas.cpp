#include <TSQR/Tsqr_Blas.hpp>
#include <TSQR/Tsqr_FortranCInterface.hpp>
#include <complex>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// mfh 28 Apr 2010
//
// C doesn't allow 'extern "C"' declarations inside a class' member
// functions, which means I have to list them all up here.

extern "C" void FortranCInterface_GLOBAL(dgemv, DGEMV) 
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

extern "C" void FortranCInterface_GLOBAL(sgemv, SGEMV)
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

extern "C" void FortranCInterface_GLOBAL(zgemv, ZGEMV)
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

extern "C" void FortranCInterface_GLOBAL(cgemv, CGEMV)
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

extern "C" void FortranCInterface_GLOBAL(dgemm, DGEMM)
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

extern "C" void FortranCInterface_GLOBAL(sgemm, SGEMM)
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

extern "C" void FortranCInterface_GLOBAL(zgemm, ZGEMM)
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

extern "C" void FortranCInterface_GLOBAL(cgemm, CGEMM)
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

extern "C" void FortranCInterface_GLOBAL(dger, DGER)
  (const int* const M,
   const int* const N,
   const double* const ALPHA,
   const double X[],
   const int* const INCX,
   const double Y[],
   const int* const INCY,
   double A[],
   const int* const LDA);

extern "C" void FortranCInterface_GLOBAL(sger, SGER)
  (const int* const M,
   const int* const N,
   const float* const ALPHA,
   const float X[],
   const int* const INCX,
   const float Y[],
   const int* const INCY,
   float A[],
   const int* const LDA);

extern "C" void FortranCInterface_GLOBAL(zgerc, ZGERC)
  (const int* const M,
   const int* const N,
   const std::complex<double>* const ALPHA,
   const std::complex<double> X[],
   const int* const INCX,
   const std::complex<double> Y[],
   const int* const INCY,
   std::complex<double> A[],
   const int* const LDA);

extern "C" void FortranCInterface_GLOBAL(cgerc, CGERC)
  (const int* const M,
   const int* const N,
   const std::complex<float>* const ALPHA,
   const std::complex<float> X[],
   const int* const INCX,
   const std::complex<float> Y[],
   const int* const INCY,
   std::complex<float> A[],
   const int* const LDA);

extern "C" void FortranCInterface_GLOBAL(dtrsm, DTRSM)
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

extern "C" void FortranCInterface_GLOBAL(strsm, STRSM)
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

extern "C" void FortranCInterface_GLOBAL(ztrsm, ZTRSM)
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

extern "C" void FortranCInterface_GLOBAL(ctrsm, CTRSM)
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
  BLAS<int, double>::GEMV (const char* const trans, 
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
    FortranCInterface_GLOBAL(dgemv, DGEMV) (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, float>::GEMV (const char* const trans, 
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
    FortranCInterface_GLOBAL(sgemv, SGEMV) (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::GEMV (const char* const trans, 
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
    FortranCInterface_GLOBAL(zgemv, ZGEMV) (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::GEMV (const char* const trans, 
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
    FortranCInterface_GLOBAL(cgemv, CGEMV) (trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  template<>
  void
  BLAS<int, double>::GEMM (const char* const transa,
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
    FortranCInterface_GLOBAL(dgemm, DGEMM) (transa, transb, &m, &n, &k, &alpha,
					    A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, float>::GEMM (const char* const transa,
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
    FortranCInterface_GLOBAL(sgemm, SGEMM) (transa, transb, &m, &n, &k, &alpha,
					    A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::GEMM (const char* const transa,
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
    FortranCInterface_GLOBAL(zgemm, ZGEMM) (transa, transb, &m, &n, &k, &alpha,
					    A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::GEMM (const char* const transa,
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
    FortranCInterface_GLOBAL(cgemm, CGEMM) (transa, transb, &m, &n, &k, &alpha,
					    A, &lda, B, &ldb, &beta, C, &ldc);
  }

  template<>
  void
  BLAS<int, double>::GER (const int m,
			  const int n,
			  const double alpha,
			  const double x[],
			  const int incx,
			  const double y[],
			  const int incy,
			  double A[],
			  const int lda)
  {
    FortranCInterface_GLOBAL(dger, DGER) (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, float>::GER (const int m,
			 const int n,
			 const float alpha,
			 const float x[],
			 const int incx,
			 const float y[],
			 const int incy,
			 float A[],
			 const int lda)
  {
    FortranCInterface_GLOBAL(sger, SGER) (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::GER (const int m,
					 const int n,
					 const std::complex<double> alpha,
					 const std::complex<double> x[],
					 const int incx,
					 const std::complex<double> y[],
					 const int incy,
					 std::complex<double> A[],
					 const int lda)
  {
    FortranCInterface_GLOBAL(zgerc, ZGERC) (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::GER (const int m,
					const int n,
					const std::complex<float> alpha,
					const std::complex<float> x[],
					const int incx,
					const std::complex<float> y[],
					const int incy,
					std::complex<float> A[],
					const int lda)
  {
    FortranCInterface_GLOBAL(cgerc, CGERC) (&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  template<>
  void
  BLAS<int, double >::TRSM (const char* const side,
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
    FortranCInterface_GLOBAL(dtrsm, DTRSM) (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

  template<>
  void
  BLAS<int, float >::TRSM (const char* const side,
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
    FortranCInterface_GLOBAL(strsm, STRSM) (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

  template<>
  void
  BLAS<int, std::complex<double> >::TRSM (const char* const side,
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
    FortranCInterface_GLOBAL(ztrsm, ZTRSM) (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

  template<>
  void
  BLAS<int, std::complex<float> >::TRSM (const char* const side,
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
    FortranCInterface_GLOBAL(ctrsm, CTRSM) (side, uplo, transa, diag, &m, &n, &alpha, A, &lda, B, &ldb);
  }

} // namespace TSQR
