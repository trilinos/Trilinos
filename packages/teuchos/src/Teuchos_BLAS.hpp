/*Paul
// 27-May-2002 General cleanup. Checked for newNamingConvention (nothing changed).
06-August-2002 Changed to images (nothing changed).
*/

// Kris
// 06.16.03 -- Start over from scratch
// 06.16.03 -- Initial templatization (Tpetra_BLAS.cpp is no longer needed)
// 06.18.03 -- Changed xxxxx_() function calls to XXXXX_F77()
//          -- Added warning messages for default calls
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_BLAS_HPP_
#define _TEUCHOS_BLAS_HPP_

#include "Teuchos_BLAS_wrappers.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Teuchos
{
//! Teuchos::BLAS: The Templated Petra BLAS Class.
/*! The Teuchos::BLAS class provides functionality similar to the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vectoer multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Teuchos_BLAS class provides C++ bindings for the BLAS
    kernels in order to insulate the rest of Petra from the details of 
    C++ to Fortran translation.

    In addition to giving access the standard BLAS functionality.
    Teuchos::BLAS also provide functionality for any <scalarType> class that
    defines the +, - * and / operators.

    Teuchos::BLAS is a single memory image interface only.  This is appropriate 
    since the standard BLAS are only specified for serial execution 
    (or shared memory parallel).
*/

  // TRANS  = ` No transpose', ` Transpose', ` Conjugate transpose' ( X, X T, XC )
  // UPLO  = ` Upper triangular', ` Lower triangular'
  // DIAG  = ` Non-unit triangular', ` Unit triangular'
  // SIDE  = ` Left', ` Right' (A or op(A) on the left, or A or op(A) on the right)

  template<typename OrdinalType, typename ScalarType>
  class BLAS
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    ScalarType ASUM(int, ScalarType*, int);
    void AXPY(int, ScalarType, ScalarType*, int, ScalarType*, int);
    void COPY(int, ScalarType*, int, ScalarType*, int);
    ScalarType DOT(int, ScalarType*, int, ScalarType*, int);
    ScalarType NRM2(int, ScalarType*, int);
    void SCAL(int, ScalarType, ScalarType*, int);
    int IAMAX(int, ScalarType*, int);
    void GEMV(char, int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType, ScalarType*, int);
    void TRMV(char, char, char, int, ScalarType*, int, ScalarType*, int);
    void GER(int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType*, int);
    void GEMM(char, char, int, int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType, ScalarType*, int);
    void SYMM(char, char, int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType beta, ScalarType*, int);
    void TRMM(char, char, char, char, int, int, ScalarType, ScalarType*, int, ScalarType*, int);
    void TRSM(char, char, char, char, int*, int*, ScalarType*, ScalarType*, int*, ScalarType*, int*);
    void XERBLA(char, int);
  };

  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::ASUM(int n, ScalarType* x, int incx)
  {
    std::cout << "Warning: default BLAS::ASUM() not yet implemented" << std::endl;
    ScalarType dummy;
    return dummy;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::AXPY(int n, ScalarType alpha, ScalarType* x, int incx, ScalarType* y, int incy)
  {
    std::cout << "Warning: default BLAS::AXPY() still experimental" << std::endl;
    if((incx == incy) && (incx > 0))
      {
	int i;
	for(i = 0; i < n; i += incx)
	  {
	    y[i] += alpha * x[i];
	  }
      }
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::COPY(int n, ScalarType* x, int incx, ScalarType* y, int incy)
  {
    std::cout << "Warning: default BLAS::COPY() still experimental" << std::endl;
    if((incx == incy) && (incx > 0))
      {
	int i;
	for(i = 0; i < n; i += incx)
	  {
	    y[i] = x[i];
	  }
      }
  }
  
  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::DOT(int n, ScalarType* x, int incx, ScalarType* y, int incy)
  {
    std::cout << "Warning: default BLAS::DOT() still experimental" << std::endl;
    ScalarType result = ScalarTraits<ScalarType>::zero();
    if((incx == incy) && (incx > 0))
      {
	int i;
	for(i = 0; i < n; i+= incx)
	  {
	    result += x[i] * y[i];
	  }
      }
    return result;
  }
  
  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::NRM2(int n, ScalarType* x, int incx)
  {
    std::cout << "Warning: default BLAS::NRM2() still experimental" << std::endl;
    ScalarType result = ScalarTraits<ScalarType>::zero();
    int i;
    for(i = 0; i < n; i += incx)
      {
	result += x[i] * x[i];
      }
    // what if sqrt() is not defined for ScalarType?
    result = sqrt(result);
    return result;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SCAL(int n, ScalarType alpha, ScalarType* x, int incx)
  {
    std::cout << "Warning: default BLAS::SCAL() still experimental" << std::endl;
    int i;
    for(i = 0; i < n; i += incx)
      {
	x[i] *= alpha;
      }
  }
  
  template<typename OrdinalType, typename ScalarType>
  int BLAS<OrdinalType, ScalarType>::IAMAX(int n, ScalarType* x, int incx)
  {
    std::cout << "Warning: default BLAS::IAMAX() not yet implemented" << std::endl;
    ScalarType dummy;
    return dummy;
  }

// TRANS  = ` No transpose', ` Transpose', ` Conjugate transpose' ( X, X T, XC )
// UPLO  = ` Upper triangular', ` Lower triangular'
// DIAG  = ` Non-unit triangular', ` Unit triangular'
// SIDE  = ` Left', ` Right' (A or op(A) on the left, or A or op(A) on the right)

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GEMV(char trans, int m, int n, ScalarType alpha, ScalarType* A, int lda, ScalarType* x, int incx, ScalarType beta, ScalarType* y, int incy)
  {
    // y <- aAx + by (A = A, A', A^H)
    std::cout << "Warning: default BLAS::GEMV() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMV(char uplo, char trans, char diag, int n, ScalarType* a, int lda, ScalarType* x, int incx)
  {
    // x <- Ax (A = A, A', A^H)
    std::cout << "Warning: default BLAS::TRMV() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GER(int m, int n, ScalarType alpha, ScalarType* x, int incx, ScalarType* y, int incy, ScalarType* A, int lda)
  {
    // A <- axy' + A
    std::cout << "Warning: default BLAS::GER() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GEMM(char transa, char transb, int m, int n, int k, ScalarType alpha, ScalarType* a, int lda, ScalarType* b, int ldb, ScalarType beta, ScalarType* c, int ldc)
  {
    std::cout << "Warning: default BLAS::GEMM() still experimental" << std::endl;
    int incra, incca, incrb, inccb;
    if (transa=='N')
      {
	incra = lda; 
	incca = 1;
      }
    else
      {
	incca = lda; 
	incra = 1;
      }
    if (transb=='N')
      {
	incrb = ldb; 
	inccb = 1;
      }
    else
      {
	inccb = ldb; 
	incrb = 1;
      }
    ScalarType* curC = c;
    ScalarType* curB = b;
    ScalarType* curA = a;
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one =  ScalarTraits<ScalarType>::one();
    if (beta==zero) 
      for (int i=0; i<m; i++) 
	c[i] = zero;
    else 
      for (int i=0; i<m; i++) 
	c[i] *= beta;
    
    for (int i=0; i<m; i++)
      {
	for (int j=0; j<n; j++)
	  {
	    ScalarType* tmpB = curB;
	    ScalarType* tmpA = curA;
	    ScalarType  tmpC = zero;
	    for (int l=0; l<k; l++)
	      {
		tmpC += (*tmpA) * (*tmpB); // Note: This is not optimal.  work on it later
		tmpA += incra;
		tmpB += inccb;
	      }
	    *curC += alpha * tmpC;
	    curC += ldc;
	    curB += incrb;
	  }
	curA += incca;
	curB = b;
	curC = c+i+1;
      }
    return;
  }
  
// TRANS  = ` No transpose', ` Transpose', ` Conjugate transpose' ( X, X T, XC )
// UPLO  = ` Upper triangular', ` Lower triangular'
// DIAG  = ` Non-unit triangular', ` Unit triangular'
// SIDE  = ` Left', ` Right' (A or op(A) on the left, or A or op(A) on the right)

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SYMM(char side, char uplo, int m, int n, ScalarType alpha, ScalarType* a, int lda, ScalarType* b, int ldb, ScalarType beta, ScalarType* c, int ldc)
  {
    // C <- aAB + bC || C <- aBA + bC | A = A'
    // UPLO refers to characteristic of B?
    std::cout << "Warning: default BLAS::SYMM() not yet implemented" << std::endl;
    if(side == 'L')
      {
	if(uplo == 'U')
	  {

	  }
	else if(uplo == 'L')
	  {

	  }
	else
	  {
	    // bad 'UPLO' parameter
	  }
      }
    else if(side == 'R')
      {
	if(uplo == 'U')
	  {

	  }
	else if(uplo == 'L')
	  {

	  }
	else
	  {
	    // bad 'UPLO' parameter
	  }
      }
    else
      {
	// bad 'SIDE' parameter
      }
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMM(char side, char uplo, char transa, char diag, int m, int n, ScalarType alpha, ScalarType* a, int lda, ScalarType* b, int ldb)
  {
    // B <- aAB || B <- aBA | A = A || A' || A^H
    std::cout << "Warning: default BLAS::TRMM() not yet implemented" << std::endl; 
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRSM(char side, char uplo, char transa, char diag, int* m, int* n, ScalarType* alpha, ScalarType* a, int* lda, ScalarType* b, int* ldb)
  {
    // B <- aA^-1B, B <- aBA^-1 | A = A || A' || A^H
    std::cout << "Warning: default BLAS::TRSM() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::XERBLA(char xerbla_arg, int info)
  {
    std::cout << "Warning: default BLAS::XERBLA() not yet implemented" << std::endl;
  }

  template<typename OrdinalType>
  class BLAS<OrdinalType, float>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    float ASUM(int, float*, int);
    void AXPY(int, float, float*, int, float*, int);
    void COPY(int, float*, int, float*, int);
    float DOT(int, float*, int, float*, int);
    float NRM2(int, float*, int);
    void SCAL(int, float, float*, int);
    int IAMAX(int, float*, int);
    void GEMV(char, int, int, float, float*, int, float*, int, float, float*, int);
    void TRMV(char, char, char, int, float*, int, float*, int);
    void GER(int, int, float, float*, int, float*, int, float*, int);
    void GEMM(char, char, int, int, int, float, float*, int, float*, int, float, float*, int);
    void SYMM(char, char, int, int, float, float*, int, float*, int, float, float*, int);
    void TRMM(char, char, char, char, int, int, float, float*, int, float*, int);
    void TRSM(char, char, char, char, int*, int*, float*, float*, int*, float*, int*);
    void XERBLA(char, int);
  };

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::ASUM(int n, float* x, int incx)
  {
    return SASUM_F77(&n, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::AXPY(int n, float alpha, float* x, int incx, float* y, int incy)
  {
    SAXPY_F77(&n, &alpha, x, &incx, y, &incy);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::COPY(int n, float* x, int incx, float* y, int incy)
  {
    SCOPY_F77(&n, x, &incx, y, &incy);
  }
  
  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::DOT(int n, float* x, int incx, float* y, int incy)
  {
    return SDOT_F77(&n, x, &incx, y, &incy);
  }
  
  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::NRM2(int n, float* x, int incx)
  {
    return SNRM2_F77(&n, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SCAL(int n, float alpha, float* x, int incx)
  {
    SSCAL_F77(&n, &alpha, x, &intx);
  }
  
  template<typename OrdinalType>
  int BLAS<OrdinalType, float>::IAMAX(int n, float* x, int incx)
  {
    return ISAMAX_F77(&n, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMV(char trans, int m, int n, float alpha, float* A, int lda, float* x, int incx, float beta, float* y, int incy)
  {
    SGEMV_F77(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMV(char uplo, char trans, char diag, int n, float* a, int lda, float* x, int incx)
  {
    STRMV_F77(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GER(int m, int n, float alpha, float* x, int incx, float* y, int incy, float* a, int lda)
  {
    SGER_F77(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMM(char transa, char transb, int m, int n, int k, float alpha, float* a, int lda, float* b, int ldb, float beta, float* c, int ldc)
  {
    SGEMM_F77(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SYMM(char side, char uplo, int m, int n, float alpha, float* a, int lda, float *b, int ldb, float beta, float *c, int ldc)
  {
    SSYMM_F77(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMM(char side, char uplo, char transa, char diag, int m, int n, float alpha, float* a, int lda, float* b, int ldb)
  {
    STRMM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRSM(char side, char uplo, char transa, char diag, int* m, int* n, float* alpha, float* a, int* lda, float* b, int* ldb)
  {
    STRSM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::XERBLA(char xerbla_arg, int info)
  {
    XERBLA_F77(&xerbla_arg, &info);
  }

  template<typename OrdinalType>
  class BLAS<OrdinalType, double>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    double ASUM(int, double* x, int);
    void AXPY(int, double, double*, int, double*, int);
    void COPY(int, double*, int, double*, int);
    double DOT(int, double* x, int, double* y, int);
    double NRM2(int, double* x, int);
    void SCAL(int, double, double*, int);
    int IAMAX(int, double*, int);
    void GEMV(char, int, int, double, double*, int, double*, int, double, double*, int);
    void TRMV(char, char, char, int, double*, int, double*, int);
    void GER(int, int, double, double*, int, double*, int, double*, int);
    void GEMM(char, char, int, int, int, double, double*, int, double*, int, double, double*, int);
    void SYMM(char, char, int, int, double, double*, int, double*, int, double, double*, int);
    void TRMM(char, char, char, char, int, int, double, double*, int, double*, int);
    void TRSM(char, char, char, char, int*, int*, double*, double*, int*, double*, int*);
    void XERBLA(char, int);
  };
    
  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::ASUM(int n, double* x, int incx)
  {
    return DASUM_F77(&n, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::AXPY(int n, double alpha, double* x, int incx, double* y, int incy)
  {
    DAXPY_F77(&n, &alpha, x, &incx, y, &incy);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::COPY(int n, double* x, int incx, double* y, int incy)
  {
    DCOPY_F77(&n, x, &incx, y, &incy);
  }
  
  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::DOT(int n, double* x, int incx, double* y, int incy)
  {
    return DDOT_F77(&n, x, &incx, y, &incy);
  }
  
  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::NRM2(int n, double* x, int incx)
  {
    return DNRM2_F77(&n, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SCAL(int n, double alpha, double* x, int incx)
  {
    DSCAL_F77(&n, &alpha, x, &incx);
  }
  
  template<typename OrdinalType>
  int BLAS<OrdinalType, double>::IAMAX(int n, double* x, int incx)
  {
    return IDAMAX_F77(&n, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMV(char trans, int m, int n, double alpha, double* A, int lda, double* x, int incx, double beta, double* y, int incy)
  {
    DGEMV_F77(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMV(char uplo, char trans, char diag, int n, double* a, int lda, double* x, int incx)
  {
    DTRMV_F77(&uplo, trans, &diag, &n, a, &lda, x, &incx);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GER(int m, int n, double alpha, double* x, int incx, double* y, int incy, double* a, int lda)
  {
    DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMM(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc)
  {
    DGEMM_F77(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SYMM(char side, char uplo, int m, int n, double alpha, double* a, int lda, double *b, int ldb, double beta, double *c, int ldc)
  {
    DSYMM_F77(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b, int ldb)
  {
    DTRMM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRSM(char side, char uplo, char transa, char diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb)
  {
    DTRSM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
  }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::XERBLA(char xerbla_arg, int info)
  {
    XERBLA_F77(&xerbla_arg, &info);
  }
  
} // end of namespace Teuchos

#endif // end of _TEUCHOS_BLAS_HPP_
