
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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

static int epetra_blas_one = 1;

//=============================================================================
float Epetra_BLAS::ASUM(int N, float * X) const {
 return(SASUM_F77(&N, X, &epetra_blas_one));
}
//=============================================================================
double Epetra_BLAS::ASUM(int N, double * X) const {
  return(DASUM_F77(&N, X, &epetra_blas_one));
}
//=============================================================================
float Epetra_BLAS::DOT(int N, float * X, float * Y) const {
  return(SDOT_F77(&N, X, &epetra_blas_one, Y, &epetra_blas_one));
}
//=============================================================================
double Epetra_BLAS::DOT(int N, double * X, double * Y) const {
  return(DDOT_F77(&N, X, &epetra_blas_one, Y, &epetra_blas_one));
}
//=============================================================================
float Epetra_BLAS::NRM2(int N, float * X) const {
  return(SNRM2_F77(&N, X, &epetra_blas_one));
}
//=============================================================================
double Epetra_BLAS::NRM2(int N, double * X) const {
  return(DNRM2_F77(&N, X, &epetra_blas_one));
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, float ALPHA, float * X) const {
  SSCAL_F77(&N, &ALPHA, X, &epetra_blas_one);
  return;
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, double ALPHA, double * X) const {
  DSCAL_F77(&N, &ALPHA, X, &epetra_blas_one);
  return;
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, float * X) const {
  return(ISAMAX_F77(&N, X, &epetra_blas_one)-1);// Note that we return base zero result
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, double * X) const {
  return(IDAMAX_F77(&N, X, &epetra_blas_one)-1);// Note that we return base zero result
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, float ALPHA, float * X, float * Y) const {
  SAXPY_F77(&N, &ALPHA, X, &epetra_blas_one, Y, &epetra_blas_one);
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, double ALPHA, double * X, double * Y) const {
  DAXPY_F77(&N, &ALPHA, X, &epetra_blas_one, Y, &epetra_blas_one);
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      float ALPHA, float * A, int LDA, float * X,
		      float BETA, float * Y) const {
  SGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &epetra_blas_one, &BETA, Y, &epetra_blas_one);
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      double ALPHA, double * A, int LDA, double * X,
		      double BETA, double * Y) const {
  DGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &epetra_blas_one, &BETA, Y, &epetra_blas_one);
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
