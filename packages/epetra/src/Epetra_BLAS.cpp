
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


#include "Epetra_BLAS.h"
#include "Epetra_BLAS_wrappers.h"


//=============================================================================
float Epetra_BLAS::ASUM(int N, float * X) const {
  int one = 1;
 return(sasum_(&N, X, &one));
}
//=============================================================================
double Epetra_BLAS::ASUM(int N, double * X) const {
  int one = 1;
  return(dasum_(&N, X, &one));
}
//=============================================================================
float Epetra_BLAS::DOT(int N, float * X, float * Y) const {
  int one = 1;
  return(sdot_(&N, X, &one, Y, &one));
}
//=============================================================================
double Epetra_BLAS::DOT(int N, double * X, double * Y) const {
  int one = 1;
  return(ddot_(&N, X, &one, Y, &one));
}
//=============================================================================
float Epetra_BLAS::NRM2(int N, float * X) const {
  int one = 1;
  return(snrm2_(&N, X, &one));
}
//=============================================================================
double Epetra_BLAS::NRM2(int N, double * X) const {
  int one = 1;
  return(dnrm2_(&N, X, &one));
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, float ALPHA, float * X) const {
  int one = 1;
  sscal_(&N, &ALPHA, X, &one);
  return;
}
//=============================================================================
void Epetra_BLAS::SCAL(int N, double ALPHA, double * X) const {
  int one = 1;
  dscal_(&N, &ALPHA, X, &one);
  return;
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, float * X) const {
  int one = 1;
  return(isamax_(&N, X, &one)-1);// Note that we return base zero result
}
//=============================================================================
int Epetra_BLAS::IAMAX(int N, double * X) const {
  int one = 1;
  return(idamax_(&N, X, &one)-1);// Note that we return base zero result
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, float ALPHA, float * X, float * Y) const {
  int one = 1;
  saxpy_(&N, &ALPHA, X, &one, Y, &one);
}
//=============================================================================
void Epetra_BLAS::AXPY(int N, double ALPHA, double * X, double * Y) const {
  int one = 1;
  daxpy_(&N, &ALPHA, X, &one, Y, &one);
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      float ALPHA, float * A, int LDA, float * X,
		      float BETA, float * Y) const {
#if defined(INTEL_CXML)
  unsigned int one = 1;
  sgemv_(&TRANS, one, &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
#else
  int one = 1;
  sgemv_(&TRANS, &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
#endif
}
//=============================================================================
void Epetra_BLAS::GEMV(char TRANS, int M, int N,
		      double ALPHA, double * A, int LDA, double * X,
		      double BETA, double * Y) const {
  
#if defined(INTEL_CXML)
  unsigned int one = 1;
  dgemv_(&TRANS, one, &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
#else
  int one = 1;
  dgemv_(&TRANS, &M, &N, &ALPHA,
	 A, &LDA, X, &one, &BETA, Y, &one);
#endif
}

//=============================================================================
void Epetra_BLAS::GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const {

#if defined(INTEL_CXML)
	unsigned int one = 1;
	sgemm_(&TRANSA, one, &TRANSB, one, &M, &N, &K, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#else
	sgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#endif
}

//=============================================================================
void Epetra_BLAS::GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const {
#if defined(INTEL_CXML)
	unsigned int one = 1;
    dgemm_(&TRANSA, one, &TRANSB, one, &M, &N, &K, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#else
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#endif
}
//=============================================================================
void Epetra_BLAS::SYMM(char SIDE, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB, float BETA, float * C, int LDC) const {

#if defined(INTEL_CXML)
	unsigned int one = 1;
	ssymm_(&SIDE, one, &M, &N, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#else
	ssymm_(&SIDE, &M, &N, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#endif
}

//=============================================================================
void Epetra_BLAS::SYMM(char SIDE, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB, double BETA, double * C, int LDC) const {
#if defined(INTEL_CXML)
	unsigned int one = 1;
    dsymm_(&SIDE, one, &M, &N, &K, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#else
    dsymm_(&SIDE, &M, &N, &ALPHA,
	   A, &LDA, B, &LDB, &BETA, C, &LDC);
#endif
}
//=============================================================================
void Epetra_BLAS::TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    float ALPHA, float * A, int LDA, float * B,
	    int LDB) const {
#if defined(INTEL_CXML)
	unsigned int one = 1;
	    strmm_(&SIDE, one, &UPLO, one, &TRANSA, one, &DIAG, one, &M, &N, &ALPHA, 
			A, &LDA, B, &LDB);
#else
    strmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, A, &LDA, B, &LDB);
#endif
}
//=============================================================================
void Epetra_BLAS::TRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N,
	    double ALPHA, double * A, int LDA, double * B,
	    int LDB) const {
#if defined(INTEL_CXML)
	unsigned int one = 1;
    dtrmm_(&SIDE, one, &UPLO, one, &TRANSA, one, &DIAG, one, &M, &N, &ALPHA, A, &LDA, B, &LDB);
#else
    dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, A, &LDA, B, &LDB);
#endif
}
