
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

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
float Epetra_BLAS::ASUM(int N, float * X) const {
  int one = 1;
 return(SASUM_F77(&N, X, &one));
}
//=============================================================================
double Epetra_BLAS::ASUM(int N, double * X) const {
  int one = 1;
  return(DASUM_F77(&N, X, &one));
}
//=============================================================================
float Epetra_BLAS::DOT(int N, float * X, float * Y) const {
  int one = 1;
  return(SDOT_F77(&N, X, &one, Y, &one));
}
//=============================================================================
double Epetra_BLAS::DOT(int N, double * X, double * Y) const {
  int one = 1;
  return(DDOT_F77(&N, X, &one, Y, &one));
}
//=============================================================================
float Epetra_BLAS::NRM2(int N, float * X) const {
  int one = 1;
  return(SNRM2_F77(&N, X, &one));
}
//=============================================================================
double Epetra_BLAS::NRM2(int N, double * X) const {
  int one = 1;
  return(DNRM2_F77(&N, X, &one));
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, float ALPHA, float * X) const {
  int one = 1;
  SSCAL_F77(&N, &ALPHA, X, &one);
  return;
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, double ALPHA, double * X) const {
  int one = 1;
  DSCAL_F77(&N, &ALPHA, X, &one);
  return;
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, float * X) const {
  int one = 1;
  return(ISAMAX_F77(&N, X, &one)-1);// Note that we return base zero result
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, double * X) const {
  int one = 1;
  return(IDAMAX_F77(&N, X, &one)-1);// Note that we return base zero result
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, float ALPHA, float * X, float * Y) const {
  int one = 1;
  SAXPY_F77(&N, &ALPHA, X, &one, Y, &one);
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, double ALPHA, double * X, double * Y) const {
  int one = 1;
  DAXPY_F77(&N, &ALPHA, X, &one, Y, &one);
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      float ALPHA, float * A, int LDA, float * X,
		      float BETA, float * Y) const {
  int one = 1;
  SGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      double ALPHA, double * A, int LDA, double * X,
		      double BETA, double * Y) const {
  int one = 1;
  DGEMV_F77(CHAR_MACRO(TRANS), &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
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
