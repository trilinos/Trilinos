/*Paul
27-May-2002 General cleanup. Changed method names to fit namingConvention (nothing changed).
06-August-2002 Changed to images (nothing changed).
*/

#if defined (INTEL_CXML)
#define UPLO_MACRO &UPLO, one
#define TRANS_MACRO &TRANS, one
#define FACT_MACRO &FACT, one
#define EQUED_MACRO &EQUED, one
#define NORM_MACRO &NORM, one
#define JOB_MACRO &JOB, one
#define COMPZ_MACRO &COMPZ, one
#define SIDE_MACRO &SIDE, one
#define HOWMNY_MACRO &HOWMNY, one
#define COMPQ_MACRO &COMPQ, one
#define CMACH_MACRO &CMACH, one
#else
#define UPLO_MACRO &UPLO
#define TRANS_MACRO &TRANS
#define FACT_MACRO &FACT
#define EQUED_MACRO &EQUED
#define NORM_MACRO &NORM
#define JOB_MACRO &JOB
#define COMPZ_MACRO &COMPZ
#define SIDE_MACRO &SIDE
#define HOWMNY_MACRO &HOWMNY
#define COMPQ_MACRO &COMPQ
#define CMACH_MACRO &CMACH
#endif

#include "Tpetra_LAPACK_wrappers.h"

namespace Tpetra 
{
// Symmetric positive definite linear systems

//=============================================================================
void LAPACK::POTRF( char UPLO, int N, float * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPOTRF_F77(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void LAPACK::POTRF( char UPLO, int N, double * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPOTRF_F77(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void LAPACK::POTRS( char UPLO, int N, int NRHS, float * A, int LDA, 
			  float * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPOTRS_F77(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::POTRS( char UPLO, int N, int NRHS, double * A, int LDA, 
			  double * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPOTRS_F77(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::POTRI( char UPLO, int N, float * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPOTRI_F77(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void LAPACK::POTRI( char UPLO, int N, double * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPOTRI_F77(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void LAPACK::POCON( char UPLO, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPOCON_F77(UPLO_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::POCON( char UPLO, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPOCON_F77(UPLO_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::POSV( char UPLO, int N, int NRHS, float * A, int LDA, 
			  float * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPOSV_F77(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::POSV( char UPLO, int N, int NRHS, double * A, int LDA, 
			  double * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPOSV_F77(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::POEQU(int N, float * A, int LDA, float * S, float * SCOND, 
			 float * AMAX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void LAPACK::POEQU(int N, double * A, int LDA, double * S, double * SCOND,
			double * AMAX, int * INFO) const {		 
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void LAPACK::PORFS(char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     float * B, int LDB, float * X, int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPORFS_F77(UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::PORFS(char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPORFS_F77( UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF,B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::POSVX(char FACT, char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     char EQUED, float * S, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SPOSVX_F77(FACT_MACRO, UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, EQUED_MACRO, S, B, &LDB, X, &LDX, 
	  RCOND, FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::POSVX(char FACT, char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     char EQUED, double * S, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DPOSVX_F77(FACT_MACRO, UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, EQUED_MACRO, S, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}

// General linear systems
//=============================================================================
void LAPACK::GELS( char TRANS, int m, int n, int numrhs,
						double* a, int lda, double* b, int ldb, double* work, int lwork, 
						int info) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
	DGELS_F77 (TRANS_MACRO, &m, &n, &numrhs, a, &lda, b, &ldb, work, &lwork, &info );
				
}
//=============================================================================
void LAPACK::GETRF( int M, int N, float * A, int LDA, int * IPIV, int * INFO) const {
  SGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void LAPACK::GETRF( int M, int N, double * A, int LDA, int * IPIV, int * INFO) const {
  DGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void LAPACK::GETRS( char TRANS, int N, int NRHS, float * A, int LDA, 
			  int * IPIV, float * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SGETRS_F77(TRANS_MACRO, &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::GETRS( char TRANS, int N, int NRHS, double * A, int LDA, 
			  int * IPIV, double * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DGETRS_F77(TRANS_MACRO, &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::GETRI( int N, float * A, int LDA, int * IPIV, 
			  float * WORK, int * LWORK, int * INFO) const {
  SGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void LAPACK::GETRI( int N, double * A, int LDA, int * IPIV, 
			  double * WORK, int * LWORK, int * INFO) const {
  DGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void LAPACK::GECON( char NORM, int N, float  * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SGECON_F77(NORM_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::GECON( char NORM, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DGECON_F77(NORM_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::GESV( int N, int NRHS, float * A, int LDA, int * IPIV, 
			  float * X, int LDX, int * INFO) const {
  SGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::GESV( int N, int NRHS, double * A, int LDA, int * IPIV, 
			  double * X, int LDX, int * INFO) const {
  DGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void LAPACK::GEEQU(int M, int N, float * A, int LDA, float * R, float * C, 
			 float * ROWCND, float * COLCND, float * AMAX, int * INFO) const {
  SGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void LAPACK::GEEQU(int M, int N, double * A, int LDA, double * R, double * C,  
			 double * ROWCND, double * COLCND, double * AMAX, int * INFO) const {
  DGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void LAPACK::GERFS(char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	      int * IPIV, float * B, int LDB, float * X, int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SGERFS_F77(TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::GERFS(char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     int * IPIV, double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DGERFS_F77( TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::GESVX(char FACT, char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
			 int * IPIV, char EQUED, float * R, float * C, float * B, int LDB, float * X, int LDX, float * RCOND, 
			 float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SGESVX_F77(FACT_MACRO, TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, EQUED_MACRO, R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void LAPACK::GESVX(char FACT, char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
			 int * IPIV, char EQUED, double * R, double * C, double * B, int LDB, double * X, int LDX, double * RCOND, 
			 double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
 #if defined (INTEL_CXML)
	unsigned int one=1;
#endif
 DGESVX_F77(FACT_MACRO, TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, EQUED_MACRO, R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}




//=============================================================================
void LAPACK::GEHRD(int N, int ILO, int IHI, float * A, int LDA, float * TAU, 
			 float * WORK, int LWORK, int * INFO) const {
  SGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::GEHRD(int N, int ILO, int IHI, double * A, int LDA, double * TAU, 
			 double * WORK, int LWORK, int * INFO) const {
  DGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, float * H, int LDH, 
			  float * WR, float * WI, float * Z, int LDZ, float * WORK, int LWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SHSEQR_F77(JOB_MACRO, COMPZ_MACRO, &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, double * H, int LDH, 
			  double * WR, double * WI, double * Z, int LDZ, double * WORK, int LWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DHSEQR_F77(JOB_MACRO, COMPZ_MACRO, &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::ORGHR( int N, int ILO, int IHI, float * A, int LDA, float * TAU, 
			  float * WORK, int LWORK, int * INFO) const {
  SORGHR_F77( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::ORGHR( int N, int ILO, int IHI, double * A, int LDA, double * TAU, 
			  double * WORK, int LWORK, int * INFO) const {
  DORGHR_F77( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, float * A, int LDA, 
			  float * TAU, float * C, int LDC, float * WORK, int LWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  SORMHR_F77(SIDE_MACRO, TRANS_MACRO, &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, double * A, int LDA, 
			  double * TAU, double * C, int LDC, double * WORK, int LWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DORMHR_F77(SIDE_MACRO, TRANS_MACRO, &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void LAPACK::TREVC( char SIDE, char HOWMNY, int * SELECT, int N, float * T, int LDT, float *VL, int LDVL,
			  float * VR, int LDVR, int MM, int * M, float * WORK, int * INFO) const {

#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  STREVC_F77(SIDE_MACRO, HOWMNY_MACRO, SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void LAPACK::TREVC( char SIDE, char HOWMNY, int * SELECT, int N, double * T, int LDT, double *VL, int LDVL,
			  double * VR, int LDVR, int MM, int * M, double * WORK, int * INFO) const {

#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  DTREVC_F77(SIDE_MACRO, HOWMNY_MACRO, SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void LAPACK::TREXC( char COMPQ, int N, float * T, int LDT, float * Q, int LDQ, int IFST, int ILST, 
			  float * WORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  STREXC_F77( COMPQ_MACRO, &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
void LAPACK::TREXC( char COMPQ, int N, double * T, int LDT, double * Q, int LDQ, int IFST, int ILST, 
			  double * WORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  DTREXC_F77( COMPQ_MACRO, &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
float LAPACK::SLAMCH( char CMACH) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  return(SLAMCH_F77( CMACH_MACRO));
}
//=============================================================================
double LAPACK::DLAMCH( char CMACH) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  return(DLAMCH_F77( CMACH_MACRO));
}
  
} // namespace Tpetra
