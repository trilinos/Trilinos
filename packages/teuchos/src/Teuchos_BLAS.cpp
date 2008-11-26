// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Teuchos_BLAS.hpp"
#include "Teuchos_BLAS_wrappers.hpp"

#ifdef TEUCHOS_BLAS_APPLE_VECLIB_ERROR
#include <vecLib/cblas.h>
#endif

const char Teuchos::ESideChar[] = {'L' , 'R' };
const char Teuchos::ETranspChar[] = {'N' , 'T' , 'C' };
const char Teuchos::EUploChar[] = {'U' , 'L' };
const char Teuchos::EDiagChar[] = {'U' , 'N' };
//const char Teuchos::EFactChar[] = {'F', 'N' };
//const char Teuchos::ENormChar[] = {'O', 'I' };
//const char Teuchos::ECompQChar[] = {'N', 'I', 'V' };
//const char Teuchos::EJobChar[] = {'E', 'V', 'B' };
//const char Teuchos::EJobSChar[] = {'E', 'S' };
//const char Teuchos::EJobVSChar[] = {'V', 'N' };
//const char Teuchos::EHowmnyChar[] = {'A', 'S' };
//const char Teuchos::ECMachChar[] = {'E', 'S', 'B', 'P', 'N', 'R', 'M', 'U', 'L', 'O' };
//const char Teuchos::ESortChar[] = {'N', 'S'};

namespace Teuchos {

#ifdef HAVE_TEUCHOS_BLASFLOAT

  // *************************** BLAS<int,float> DEFINITIONS ******************************  

  void BLAS<int, float>::ROTG(float* da, float* db, float* c, float* s) const
  { SROTG_F77(da, db, c, s ); }

  void BLAS<int, float>::ROT(const int n, float* dx, const int incx, float* dy, const int incy, float* c, float* s) const
  { SROT_F77(&n, dx, &incx, dy, &incy, c, s); }

  
  float BLAS<int, float>::ASUM(const int n, const float* x, const int incx) const
  {
    float tmp = SASUM_F77(&n, x, &incx);
    return tmp;
  }
    
  void BLAS<int, float>::AXPY(const int n, const float alpha, const float* x, const int incx, float* y, const int incy) const
  { SAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  void BLAS<int, float>::COPY(const int n, const float* x, const int incx, float* y, const int incy) const 
  { SCOPY_F77(&n, x, &incx, y, &incy); }
  
  float BLAS<int, float>::DOT(const int n, const float* x, const int incx, const float* y, const int incy) const
  { return SDOT_F77(&n, x, &incx, y, &incy); }
  
  int BLAS<int, float>::IAMAX(const int n, const float* x, const int incx) const
  { return ISAMAX_F77(&n, x, &incx); }

  float BLAS<int, float>::NRM2(const int n, const float* x, const int incx) const
  { return SNRM2_F77(&n, x, &incx); }
  
  void BLAS<int, float>::SCAL(const int n, const float alpha, float* x, const int incx) const
  { SSCAL_F77(&n, &alpha, x, &incx); }
  
  void BLAS<int, float>::GEMV(ETransp trans, const int m, const int n, const float alpha, const float* A, const int lda, const float* x, const int incx, const float beta, float* y, const int incy) const
  { SGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  void BLAS<int, float>::GER(const int m, const int n, const float alpha, const float* x, const int incx, const float* y, const int incy, float* A, const int lda) const
  { SGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, float>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const float* A, const int lda, float* x, const int incx) const
  { STRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }
  
  void BLAS<int, float>::GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const float alpha, const float* A, const int lda, const float* B, const int ldb, const float beta, float* C, const int ldc) const
  { SGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  void BLAS<int, float>::SYMM(ESide side, EUplo uplo, const int m, const int n, const float alpha, const float* A, const int lda, const float* B, const int ldb, const float beta, float* C, const int ldc) const
  { SSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  void BLAS<int, float>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const float alpha, const float* A, const int lda, float* B, const int ldb) const
  { STRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }
  
  void BLAS<int, float>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const float alpha, const float* A, const int lda, float* B, const int ldb) const
  { STRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

#endif // HAVE_TEUCHOS_BLASFLOAT

  // *************************** BLAS<int,double> DEFINITIONS ******************************  
  
  void BLAS<int, double>::ROTG(double* da, double* db, double* c, double* s) const
  { DROTG_F77(da, db, c, s); }

  void BLAS<int, double>::ROT(const int n, double* dx, const int incx, double* dy, const int incy, double* c, double* s) const
  { DROT_F77(&n, dx, &incx, dy, &incy, c, s); }

  double BLAS<int, double>::ASUM(const int n, const double* x, const int incx) const
  { return DASUM_F77(&n, x, &incx); }
  
  void BLAS<int, double>::AXPY(const int n, const double alpha, const double* x, const int incx, double* y, const int incy) const
  { DAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  void BLAS<int, double>::COPY(const int n, const double* x, const int incx, double* y, const int incy) const
  { DCOPY_F77(&n, x, &incx, y, &incy); }
  
  double BLAS<int, double>::DOT(const int n, const double* x, const int incx, const double* y, const int incy) const
  { return DDOT_F77(&n, x, &incx, y, &incy); }
  
  int BLAS<int, double>::IAMAX(const int n, const double* x, const int incx) const
  { return IDAMAX_F77(&n, x, &incx); }

  double BLAS<int, double>::NRM2(const int n, const double* x, const int incx) const
  { return DNRM2_F77(&n, x, &incx); }
  
  void BLAS<int, double>::SCAL(const int n, const double alpha, double* x, const int incx) const
  { DSCAL_F77(&n, &alpha, x, &incx); }
  
  void BLAS<int, double>::GEMV(ETransp trans, const int m, const int n, const double alpha, const double* A, const int lda, const double* x, const int incx, const double beta, double* y, const int incy) const
  { DGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  void BLAS<int, double>::GER(const int m, const int n, const double alpha, const double* x, const int incx, const double* y, const int incy, double* A, const int lda) const
  { DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, double>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const double* A, const int lda, double* x, const int incx) const
  { DTRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }
  
  void BLAS<int, double>::GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const double alpha, const double* A, const int lda, const double* B, const int ldb, const double beta, double* C, const int ldc) const
  { DGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  void BLAS<int, double>::SYMM(ESide side, EUplo uplo, const int m, const int n, const double alpha, const double* A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) const
  { DSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  void BLAS<int, double>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const double alpha, const double* A, const int lda, double* B, const int ldb) const
  { DTRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  void BLAS<int, double>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const double alpha, const double* A, const int lda, double* B, const int ldb) const
  { DTRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }
  
#ifdef HAVE_TEUCHOS_COMPLEX

#ifdef HAVE_TEUCHOS_BLASFLOAT

  // *************************** BLAS<int,std::complex<float> > DEFINITIONS ******************************  

  void BLAS<int, std::complex<float> >::ROTG(std::complex<float>* da, std::complex<float>* db, float* c, std::complex<float>* s) const
  { CROTG_F77(da, db, c, s ); }

  void BLAS<int, std::complex<float> >::ROT(const int n, std::complex<float>* dx, const int incx, std::complex<float>* dy, const int incy, float* c, std::complex<float>* s) const
  { CROT_F77(&n, dx, &incx, dy, &incy, c, s); }

  float BLAS<int, std::complex<float> >::ASUM(const int n, const std::complex<float>* x, const int incx) const
  { return CASUM_F77(&n, x, &incx); }
  
  void BLAS<int, std::complex<float> >::AXPY(const int n, const std::complex<float> alpha, const std::complex<float>* x, const int incx, std::complex<float>* y, const int incy) const
  { CAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  void BLAS<int, std::complex<float> >::COPY(const int n, const std::complex<float>* x, const int incx, std::complex<float>* y, const int incy) const
  { CCOPY_F77(&n, x, &incx, y, &incy); }
  
  std::complex<float> BLAS<int, std::complex<float> >::DOT(const int n, const std::complex<float>* x, const int incx, const std::complex<float>* y, const int incy) const
  { 
#if defined(TEUCHOS_BLAS_APPLE_VECLIB_ERROR)
    std::complex<float> z;
    cblas_cdotc_sub(n,x,incx,y,incy,&z);
    return z;
#elif defined(HAVE_COMPLEX_BLAS_PROBLEM) && defined(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
    std::complex<float> z;
    CDOT_F77(&z, &n, x, &incx, y, &incy); 
    return z;
#else
    return CDOT_F77(&n, x, &incx, y, &incy); 
#endif
  }
  
  int BLAS<int, std::complex<float> >::IAMAX(const int n, const std::complex<float>* x, const int incx) const
  { return ICAMAX_F77(&n, x, &incx); }

  float BLAS<int, std::complex<float> >::NRM2(const int n, const std::complex<float>* x, const int incx) const
  { return CNRM2_F77(&n, x, &incx); }
  
  void BLAS<int, std::complex<float> >::SCAL(const int n, const std::complex<float> alpha, std::complex<float>* x, const int incx) const
  { CSCAL_F77(&n, &alpha, x, &incx); }
  
  void BLAS<int, std::complex<float> >::GEMV(ETransp trans, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, const std::complex<float>* x, const int incx, const std::complex<float> beta, std::complex<float>* y, const int incy) const
  { CGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  void BLAS<int, std::complex<float> >::GER(const int m, const int n, const std::complex<float> alpha, const std::complex<float>* x, const int incx, const std::complex<float>* y, const int incy, std::complex<float>* A, const int lda) const
  { CGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, std::complex<float> >::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const std::complex<float>* A, const int lda, std::complex<float>* x, const int incx) const
  { CTRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }
  
  void BLAS<int, std::complex<float> >::GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const std::complex<float> alpha, const std::complex<float>* A, const int lda, const std::complex<float>* B, const int ldb, const std::complex<float> beta, std::complex<float>* C, const int ldc) const
  { CGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); } 
 
  void BLAS<int, std::complex<float> >::SYMM(ESide side, EUplo uplo, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, const std::complex<float>* B, const int ldb, const std::complex<float> beta, std::complex<float>* C, const int ldc) const
  { CSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  void BLAS<int, std::complex<float> >::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb) const
  { CTRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }
  
  void BLAS<int, std::complex<float> >::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb) const
  { CTRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

#endif // HAVE_TEUCHOS_BLASFLOAT

  // *************************** BLAS<int,std::complex<double> > DEFINITIONS ******************************  

  void BLAS<int, std::complex<double> >::ROTG(std::complex<double>* da, std::complex<double>* db, double* c, std::complex<double>* s) const
  { ZROTG_F77(da, db, c, s); }

  void BLAS<int, std::complex<double> >::ROT(const int n, std::complex<double>* dx, const int incx, std::complex<double>* dy, const int incy, double* c, std::complex<double>* s) const
  { ZROT_F77(&n, dx, &incx, dy, &incy, c, s); }

  double BLAS<int, std::complex<double> >::ASUM(const int n, const std::complex<double>* x, const int incx) const
  { return ZASUM_F77(&n, x, &incx); }
  
  void BLAS<int, std::complex<double> >::AXPY(const int n, const std::complex<double> alpha, const std::complex<double>* x, const int incx, std::complex<double>* y, const int incy) const
  { ZAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  void BLAS<int, std::complex<double> >::COPY(const int n, const std::complex<double>* x, const int incx, std::complex<double>* y, const int incy) const
  { ZCOPY_F77(&n, x, &incx, y, &incy); }
  
  std::complex<double> BLAS<int, std::complex<double> >::DOT(const int n, const std::complex<double>* x, const int incx, const std::complex<double>* y, const int incy) const
  { 
#if defined(TEUCHOS_BLAS_APPLE_VECLIB_ERROR)
    std::complex<double> z;
    cblas_zdotc_sub(n,x,incx,y,incy,&z);
    return z;
#elif defined(HAVE_COMPLEX_BLAS_PROBLEM) && defined(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
    std::complex<double> z;
    ZDOT_F77(&z, &n, x, &incx, y, &incy); 
    return z;
#else
    return ZDOT_F77(&n, x, &incx, y, &incy); 
#endif
  }
  
  int BLAS<int, std::complex<double> >::IAMAX(const int n, const std::complex<double>* x, const int incx) const
  { return IZAMAX_F77(&n, x, &incx); }

  double BLAS<int, std::complex<double> >::NRM2(const int n, const std::complex<double>* x, const int incx) const
  { return ZNRM2_F77(&n, x, &incx); }
  
  void BLAS<int, std::complex<double> >::SCAL(const int n, const std::complex<double> alpha, std::complex<double>* x, const int incx) const
  { ZSCAL_F77(&n, &alpha, x, &incx); }
  
  void BLAS<int, std::complex<double> >::GEMV(ETransp trans, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, const std::complex<double>* x, const int incx, const std::complex<double> beta, std::complex<double>* y, const int incy) const
  { ZGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  void BLAS<int, std::complex<double> >::GER(const int m, const int n, const std::complex<double> alpha, const std::complex<double>* x, const int incx, const std::complex<double>* y, const int incy, std::complex<double>* A, const int lda) const
  { ZGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, std::complex<double> >::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const std::complex<double>* A, const int lda, std::complex<double>* x, const int incx) const
  { ZTRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }
  
  void BLAS<int, std::complex<double> >::GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const std::complex<double> alpha, const std::complex<double>* A, const int lda, const std::complex<double>* B, const int ldb, const std::complex<double> beta, std::complex<double>* C, const int ldc) const
  { ZGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  void BLAS<int, std::complex<double> >::SYMM(ESide side, EUplo uplo, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, const std::complex<double> *B, const int ldb, const std::complex<double> beta, std::complex<double> *C, const int ldc) const
  { ZSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  void BLAS<int, std::complex<double> >::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb) const
  { ZTRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  void BLAS<int, std::complex<double> >::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb) const
  { ZTRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }
  
#endif // HAVE_TEUCHOS_COMPLEX

}
