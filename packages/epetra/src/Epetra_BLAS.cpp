
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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

/* for INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
*/

#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, 1
#else
#define CHAR_MACRO(char_var) &char_var
#endif

#include "Epetra_BLAS.h"
#include "Epetra_BLAS_wrappers.h"


//=============================================================================
float Epetra_BLAS::ASUM(const int N, const float * X, const int INCX) const {
 return(SASUM_F77(&N, X, &INCX));
}
//=============================================================================
double Epetra_BLAS::ASUM(const int N, const double * X, const int INCX) const {
  return(DASUM_F77(&N, X, &INCX));
}
//=============================================================================
float Epetra_BLAS::DOT(const int N, const float * X, const float * Y, const int INCX, const int INCY) const {
  return(SDOT_F77(&N, X, &INCX, Y, &INCY));
}
//=============================================================================
double Epetra_BLAS::DOT(const int N, const double * X, const double * Y, const int INCX, const int INCY) const {
  return(DDOT_F77(&N, X, &INCX, Y, &INCY));
}
//=============================================================================
float Epetra_BLAS::NRM2(const int N, const float * X, const int INCX) const {
  return(SNRM2_F77(&N, X, &INCX));
}
//=============================================================================
double Epetra_BLAS::NRM2(const int N, const double * X, const int INCX) const {
  return(DNRM2_F77(&N, X, &INCX));
}
//=============================================================================
void Epetra_BLAS::SCAL(const int N, const float ALPHA, float * X, const int INCX) const {
  SSCAL_F77(&N, &ALPHA, X, &INCX);
  return;
}
//=============================================================================
void Epetra_BLAS::SCAL(const int N, const double ALPHA, double * X, const int INCX) const {
  DSCAL_F77(&N, &ALPHA, X, &INCX);
  return;
}
//=============================================================================
void Epetra_BLAS::COPY(const int N, const float * X, float * Y, const int INCX, const int INCY) const {
  SCOPY_F77(&N, X, &INCX, Y, &INCY);
  return;
}
//=============================================================================
void Epetra_BLAS::COPY(const int N, const double * X, double * Y, const int INCX, const int INCY) const {
  DCOPY_F77(&N, X, &INCX, Y, &INCY);
  return;
}
//=============================================================================
int Epetra_BLAS::IAMAX(const int N, const float * X, const int INCX) const {
  return(ISAMAX_F77(&N, X, &INCX)-1);// Note that we return base zero result
}
//=============================================================================
int Epetra_BLAS::IAMAX(const int N, const double * X, const int INCX) const {
  return(IDAMAX_F77(&N, X, &INCX)-1);// Note that we return base zero result
}
//=============================================================================
void Epetra_BLAS::AXPY(const int N, const float ALPHA, const float * X, float * Y, const int INCX, const int INCY) const {
  SAXPY_F77(&N, &ALPHA, X, &INCX, Y, &INCY);
}
//=============================================================================
void Epetra_BLAS::AXPY(const int N, const double ALPHA, const double * X, double * Y, const int INCX, const int INCY) const {
  DAXPY_F77(&N, &ALPHA, X, &INCX, Y, &INCY);
}
//=============================================================================
void Epetra_BLAS::GEMV(const char TRANS, const int M, const int N,
		      const float ALPHA, const float * A, const int LDA, const float * X,
		      const float BETA, float * Y, const int INCX, const int INCY) const {
  SGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &INCX, &BETA, Y, &INCY);
}
//=============================================================================
void Epetra_BLAS::GEMV(const char TRANS, const int M, const int N,
		      const double ALPHA, const double * A, const int LDA, const double * X,
		      const double BETA, double * Y, const int INCX, const int INCY) const {
  DGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &INCX, &BETA, Y, &INCY);
}

//=============================================================================
void Epetra_BLAS::GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K,
		       const float ALPHA, const float * A, const int LDA, const float * B,
		       const int LDB, const float BETA, float * C, const int LDC) const {

  SGEMM_F77(CHAR_MACRO(TRANSA), CHAR_MACRO(TRANSB), &M, &N, &K, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}

//=============================================================================
void Epetra_BLAS::GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K,
	    const double ALPHA, const double * A, const int LDA, const double * B,
	    const int LDB, const double BETA, double * C, const int LDC) const {

  DGEMM_F77(CHAR_MACRO(TRANSA), CHAR_MACRO(TRANSB), &M, &N, &K, &ALPHA,
	 A, &LDA, B, &LDB, &BETA, C, &LDC);
}
//=============================================================================
void Epetra_BLAS::SYMM(const char SIDE, const char UPLO, const int M, const int N,
	    const float ALPHA, const float * A, const int LDA, const float * B,
	    const int LDB, const float BETA, float * C, const int LDC) const {

  SSYMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), &M, &N, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}

//=============================================================================
void Epetra_BLAS::SYMM(const char SIDE, const char UPLO, const int M, const int N,
	    const double ALPHA, const double * A, const int LDA, const double * B,
	    const int LDB, const double BETA, double * C, const int LDC) const {

  DSYMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), &M, &N, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}
//=============================================================================
void Epetra_BLAS::TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, const int N,
	    const float ALPHA, const float * A, const int LDA, float * B,
	    const int LDB) const {

  STRMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), CHAR_MACRO(TRANSA), CHAR_MACRO(DIAG), 
	  &M, &N, &ALPHA, A, &LDA, B, &LDB);
}
//=============================================================================
void Epetra_BLAS::TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, const int N,
	    const double ALPHA, const double * A, const int LDA, double * B,
	    const int LDB) const {

  DTRMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), CHAR_MACRO(TRANSA), CHAR_MACRO(DIAG), 
	    &M, &N, &ALPHA, A, &LDA, B, &LDB);
}
//=============================================================================
void Epetra_BLAS::SYRK(const char UPLO, const char TRANS, const int N, const int K, const float ALPHA, const float *A,
                       const int LDA, const float BETA, float *C, const int LDC) const{
  SSYRK_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
}
//=============================================================================
void Epetra_BLAS::SYRK(const char UPLO, const char TRANS, const int N, const int K, const double ALPHA, const double *A,
                       const int LDA, const double BETA, double *C, const int LDC) const{
  DSYRK_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
}
