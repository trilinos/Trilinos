
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

#include "Epetra_LAPACK.h"
#include "Epetra_LAPACK_wrappers.h"


// Symmetric positive definite linear systems

//=============================================================================
void Epetra_LAPACK::POTRF( char UPLO, int N, float * A, int LDA, int * INFO) const {
  SPOTRF_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRF( char UPLO, int N, double * A, int LDA, int * INFO) const {
  DPOTRF_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRS( char UPLO, int N, int NRHS, float * A, int LDA, 
			  float * X, int LDX, int * INFO) const {
  SPOTRS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRS( char UPLO, int N, int NRHS, double * A, int LDA, 
			  double * X, int LDX, int * INFO) const {
  DPOTRS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRI( char UPLO, int N, float * A, int LDA, int * INFO) const {
  SPOTRI_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRI( char UPLO, int N, double * A, int LDA, int * INFO) const {
  DPOTRI_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POCON( char UPLO, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
  SPOCON_F77(CHAR_MACRO(UPLO), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POCON( char UPLO, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
  DPOCON_F77(CHAR_MACRO(UPLO), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSV( char UPLO, int N, int NRHS, float * A, int LDA, 
			  float * X, int LDX, int * INFO) const {
  SPOSV_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSV( char UPLO, int N, int NRHS, double * A, int LDA, 
			  double * X, int LDX, int * INFO) const {
  DPOSV_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POEQU(int N, float * A, int LDA, float * S, float * SCOND, 
			 float * AMAX, int * INFO) const {
  SPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POEQU(int N, double * A, int LDA, double * S, double * SCOND,
			double * AMAX, int * INFO) const {		 
  DPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::PORFS(char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     float * B, int LDB, float * X, int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SPORFS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::PORFS(char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
  DPORFS_F77( CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF,B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSVX(char FACT, char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     char EQUED, float * S, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, CHAR_MACRO(EQUED), S, B, &LDB, X, &LDX, 
	  RCOND, FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSVX(char FACT, char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     char EQUED, double * S, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
  DPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, CHAR_MACRO(EQUED), S, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}

// General linear systems
//=============================================================================
void Epetra_LAPACK::GELS( char TRANS, int m, int n, int numrhs,
						double* a, int lda, double* b, int ldb, double* work, int lwork, 
						int info) const {
	DGELS_F77 (CHAR_MACRO(TRANS), &m, &n, &numrhs, a, &lda, b, &ldb, work, &lwork, &info );
				
}
//=============================================================================
void Epetra_LAPACK::GETRF( int M, int N, float * A, int LDA, int * IPIV, int * INFO) const {
  SGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRF( int M, int N, double * A, int LDA, int * IPIV, int * INFO) const {
  DGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRS( char TRANS, int N, int NRHS, float * A, int LDA, 
			  int * IPIV, float * X, int LDX, int * INFO) const {
  SGETRS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRS( char TRANS, int N, int NRHS, double * A, int LDA, 
			  int * IPIV, double * X, int LDX, int * INFO) const {
  DGETRS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRI( int N, float * A, int LDA, int * IPIV, 
			  float * WORK, int * LWORK, int * INFO) const {
  SGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRI( int N, double * A, int LDA, int * IPIV, 
			  double * WORK, int * LWORK, int * INFO) const {
  DGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GECON( char NORM, int N, float  * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
  SGECON_F77(CHAR_MACRO(NORM), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GECON( char NORM, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
  DGECON_F77(CHAR_MACRO(NORM), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESV( int N, int NRHS, float * A, int LDA, int * IPIV, 
			  float * X, int LDX, int * INFO) const {
  SGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESV( int N, int NRHS, double * A, int LDA, int * IPIV, 
			  double * X, int LDX, int * INFO) const {
  DGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEEQU(int M, int N, float * A, int LDA, float * R, float * C, 
			 float * ROWCND, float * COLCND, float * AMAX, int * INFO) const {
  SGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEEQU(int M, int N, double * A, int LDA, double * R, double * C,  
			 double * ROWCND, double * COLCND, double * AMAX, int * INFO) const {
  DGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GERFS(char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	      int * IPIV, float * B, int LDB, float * X, int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SGERFS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GERFS(char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     int * IPIV, double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
  DGERFS_F77( CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESVX(char FACT, char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
			 int * IPIV, char EQUED, float * R, float * C, float * B, int LDB, float * X, int LDX, float * RCOND, 
			 float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, CHAR_MACRO(EQUED), R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESVX(char FACT, char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
			 int * IPIV, char EQUED, double * R, double * C, double * B, int LDB, double * X, int LDX, double * RCOND, 
			 double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
 DGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, CHAR_MACRO(EQUED), R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}




//=============================================================================
void Epetra_LAPACK::GEHRD(int N, int ILO, int IHI, float * A, int LDA, float * TAU, 
			 float * WORK, int LWORK, int * INFO) const {
  SGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEHRD(int N, int ILO, int IHI, double * A, int LDA, double * TAU, 
			 double * WORK, int LWORK, int * INFO) const {
  DGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, float * H, int LDH, 
			  float * WR, float * WI, float * Z, int LDZ, float * WORK, int LWORK, 
			  int * INFO) const {
  SHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, double * H, int LDH, 
			  double * WR, double * WI, double * Z, int LDZ, double * WORK, int LWORK, 
			  int * INFO) const {
  DHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORGHR( int N, int ILO, int IHI, float * A, int LDA, float * TAU, 
			  float * WORK, int LWORK, int * INFO) const {
  SORGHR_F77( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORGHR( int N, int ILO, int IHI, double * A, int LDA, double * TAU, 
			  double * WORK, int LWORK, int * INFO) const {
  DORGHR_F77( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, float * A, int LDA, 
			  float * TAU, float * C, int LDC, float * WORK, int LWORK, int * INFO) const {
  SORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, double * A, int LDA, 
			  double * TAU, double * C, int LDC, double * WORK, int LWORK, int * INFO) const {
  DORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::TREVC( char SIDE, char HOWMNY, int * SELECT, int N, float * T, int LDT, float *VL, int LDVL,
			  float * VR, int LDVR, int MM, int * M, float * WORK, int * INFO) const {

  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  STREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::TREVC( char SIDE, char HOWMNY, int * SELECT, int N, double * T, int LDT, double *VL, int LDVL,
			  double * VR, int LDVR, int MM, int * M, double * WORK, int * INFO) const {

  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  DTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::TREXC( char COMPQ, int N, float * T, int LDT, float * Q, int LDQ, int IFST, int ILST, 
			  float * WORK, int * INFO) const {
  STREXC_F77( CHAR_MACRO(COMPQ), &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::TREXC( char COMPQ, int N, double * T, int LDT, double * Q, int LDQ, int IFST, int ILST, 
			  double * WORK, int * INFO) const {
  DTREXC_F77( CHAR_MACRO(COMPQ), &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
float Epetra_LAPACK::SLAMCH( char CMACH) const {
  return(SLAMCH_F77( CHAR_MACRO(CMACH)));
}
//=============================================================================
double Epetra_LAPACK::DLAMCH( char CMACH) const {
	return(DLAMCH_F77( CHAR_MACRO(CMACH)));

}
  


