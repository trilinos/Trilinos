
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
float Epetra_BLAS::ASUM(int N, float * X, int INCX) const {
 return(SASUM_F77(&N, X, &INCX));
}
//=============================================================================
double Epetra_BLAS::ASUM(int N, double * X, int INCX) const {
  return(DASUM_F77(&N, X, &INCX));
}
//=============================================================================
float Epetra_BLAS::DOT(int N, float * X, float * Y, int INCX, int INCY) const {
  return(SDOT_F77(&N, X, &INCX, Y, &INCY));
}
//=============================================================================
double Epetra_BLAS::DOT(int N, double * X, double * Y, int INCX, int INCY) const {
  return(DDOT_F77(&N, X, &INCX, Y, &INCY));
}
//=============================================================================
float Epetra_BLAS::NRM2(int N, float * X, int INCX) const {
  return(SNRM2_F77(&N, X, &INCX));
}
//=============================================================================
double Epetra_BLAS::NRM2(int N, double * X, int INCX) const {
  return(DNRM2_F77(&N, X, &INCX));
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, float ALPHA, float * X, int INCX) const {
  SSCAL_F77(&N, &ALPHA, X, &INCX);
  return;
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, double ALPHA, double * X, int INCX) const {
  DSCAL_F77(&N, &ALPHA, X, &INCX);
  return;
}
//=============================================================================
void Epetra_BLAS::COPY(int N, float * X, float * Y, int INCX, int INCY) const {
  SCOPY_F77(&N, X, &INCX, Y, &INCY);
  return;
}
//=============================================================================
void Epetra_BLAS::COPY(int N, double * X, double * Y, int INCX, int INCY) const {
  DCOPY_F77(&N, X, &INCX, Y, &INCY);
  return;
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, float * X, int INCX) const {
  return(ISAMAX_F77(&N, X, &INCX)-1);// Note that we return base zero result
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, double * X, int INCX) const {
  return(IDAMAX_F77(&N, X, &INCX)-1);// Note that we return base zero result
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, float ALPHA, float * X, float * Y, int INCX, int INCY) const {
  SAXPY_F77(&N, &ALPHA, X, &INCX, Y, &INCY);
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, double ALPHA, double * X, double * Y, int INCX, int INCY) const {
  DAXPY_F77(&N, &ALPHA, X, &INCX, Y, &INCY);
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      float ALPHA, float * A, int LDA, float * X,
		      float BETA, float * Y, int INCX, int INCY) const {
  SGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &INCX, &BETA, Y, &INCY);
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      double ALPHA, double * A, int LDA, double * X,
		      double BETA, double * Y, int INCX, int INCY) const {
  DGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &INCX, &BETA, Y, &INCY);
}

//=============================================================================
void Epetra_BLAS::GEMM(char TRANSA, char TRANSB, int M, int N, int K,
		       float ALPHA, float * A, int LDA, float * B,
		       int LDB, float BETA, float * C, int LDC) const {

  SGEMM_F77(CHAR_MACRO(TRANSA), CHAR_MACRO(TRANSB), &M, &N, &K, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}

//=============================================================================
void Epetra_BLAS::GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const {

  DGEMM_F77(CHAR_MACRO(TRANSA), CHAR_MACRO(TRANSB), &M, &N, &K, &ALPHA,
	 A, &LDA, B, &LDB, &BETA, C, &LDC);
}
//=============================================================================
void Epetra_BLAS::SYMM(char SIDE, char UPLO, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const {

  SSYMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), &M, &N, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}

//=============================================================================
void Epetra_BLAS::SYMM(char SIDE, char UPLO, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const {

  DSYMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), &M, &N, &ALPHA,
         A, &LDA, B, &LDB, &BETA, C, &LDC);
}
//=============================================================================
void Epetra_BLAS::TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB) const {

  STRMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), CHAR_MACRO(TRANSA), CHAR_MACRO(DIAG), 
	  &M, &N, &ALPHA, A, &LDA, B, &LDB);
}
//=============================================================================
void Epetra_BLAS::TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB) const {

  DTRMM_F77(CHAR_MACRO(SIDE), CHAR_MACRO(UPLO), CHAR_MACRO(TRANSA), CHAR_MACRO(DIAG), 
	    &M, &N, &ALPHA, A, &LDA, B, &LDB);
}
