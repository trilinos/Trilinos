
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

#include "Petra_LAPACK.h"
#include "Petra_LAPACK_wrappers.h"


// Symmetric positive definite linear systems

//=============================================================================
void Petra_LAPACK::POTRF( char UPLO, int N, float * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  spotrf_(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void Petra_LAPACK::POTRF( char UPLO, int N, double * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dpotrf_(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void Petra_LAPACK::POTRS( char UPLO, int N, int NRHS, float * A, int LDA, 
			  float * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  spotrs_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::POTRS( char UPLO, int N, int NRHS, double * A, int LDA, 
			  double * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dpotrs_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::POTRI( char UPLO, int N, float * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  spotri_(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void Petra_LAPACK::POTRI( char UPLO, int N, double * A, int LDA, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dpotri_(UPLO_MACRO, &N, A, &LDA, INFO);
}
//=============================================================================
void Petra_LAPACK::POCON( char UPLO, int N, float * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  spocon_(UPLO_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::POCON( char UPLO, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dpocon_(UPLO_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::POSV( char UPLO, int N, int NRHS, float * A, int LDA, 
			  float * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sposv_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::POSV( char UPLO, int N, int NRHS, double * A, int LDA, 
			  double * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dposv_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::POEQU(int N, float * A, int LDA, float * S, float * SCOND, 
			 float * AMAX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  spoequ_(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void Petra_LAPACK::POEQU(int N, double * A, int LDA, double * S, double * SCOND,
			double * AMAX, int * INFO) const {		 
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dpoequ_(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void Petra_LAPACK::PORFS(char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     float * B, int LDB, float * X, int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sporfs_(UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::PORFS(char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dporfs_( UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF,B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::POSVX(char FACT, char UPLO, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	     char EQUED, float * S, float * B, int LDB, float * X, int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sposvx_(FACT_MACRO, UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, EQUED_MACRO, S, B, &LDB, X, &LDX, 
	  RCOND, FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::POSVX(char FACT, char UPLO, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     char EQUED, double * S, double * B, int LDB, double * X, int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dposvx_(FACT_MACRO, UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, EQUED_MACRO, S, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}

// General linear systems
//=============================================================================
void Petra_LAPACK::GELS( char TRANS, int m, int n, int numrhs,
						double* a, int lda, double* b, int ldb, double* work, int lwork, 
						int info) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
	dgels_ (TRANS_MACRO, &m, &n, &numrhs, a, &lda, b, &ldb, work, &lwork, &info );
				
}
//=============================================================================
void Petra_LAPACK::GETRF( int M, int N, float * A, int LDA, int * IPIV, int * INFO) const {
  sgetrf_(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void Petra_LAPACK::GETRF( int M, int N, double * A, int LDA, int * IPIV, int * INFO) const {
  dgetrf_(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void Petra_LAPACK::GETRS( char TRANS, int N, int NRHS, float * A, int LDA, 
			  int * IPIV, float * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sgetrs_(TRANS_MACRO, &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::GETRS( char TRANS, int N, int NRHS, double * A, int LDA, 
			  int * IPIV, double * X, int LDX, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dgetrs_(TRANS_MACRO, &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::GETRI( int N, float * A, int LDA, int * IPIV, 
			  float * WORK, int * LWORK, int * INFO) const {
  sgetri_(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GETRI( int N, double * A, int LDA, int * IPIV, 
			  double * WORK, int * LWORK, int * INFO) const {
  dgetri_(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GECON( char NORM, int N, float  * A, int LDA, float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sgecon_(NORM_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GECON( char NORM, int N, double * A, int LDA, double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dgecon_(NORM_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GESV( int N, int NRHS, float * A, int LDA, int * IPIV, 
			  float * X, int LDX, int * INFO) const {
  sgesv_(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::GESV( int N, int NRHS, double * A, int LDA, int * IPIV, 
			  double * X, int LDX, int * INFO) const {
  dgesv_(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Petra_LAPACK::GEEQU(int M, int N, float * A, int LDA, float * R, float * C, 
			 float * ROWCND, float * COLCND, float * AMAX, int * INFO) const {
  sgeequ_(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void Petra_LAPACK::GEEQU(int M, int N, double * A, int LDA, double * R, double * C,  
			 double * ROWCND, double * COLCND, double * AMAX, int * INFO) const {
  dgeequ_(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void Petra_LAPACK::GERFS(char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
	      int * IPIV, float * B, int LDB, float * X, int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sgerfs_(TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GERFS(char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
	     int * IPIV, double * B, int LDB, double * X, int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dgerfs_( TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GESVX(char FACT, char TRANS, int N, int NRHS, float * A, int LDA, float * AF, int LDAF, 
			 int * IPIV, char EQUED, float * R, float * C, float * B, int LDB, float * X, int LDX, float * RCOND, 
			 float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sgesvx_(FACT_MACRO, TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, EQUED_MACRO, R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GESVX(char FACT, char TRANS, int N, int NRHS, double * A, int LDA, double * AF, int LDAF, 
			 int * IPIV, char EQUED, double * R, double * C, double * B, int LDB, double * X, int LDX, double * RCOND, 
			 double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
 #if defined (INTEL_CXML)
	unsigned int one=1;
#endif
 dgesvx_(FACT_MACRO, TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, EQUED_MACRO, R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}




//=============================================================================
void Petra_LAPACK::GEHRD(int N, int ILO, int IHI, float * A, int LDA, float * TAU, 
			 float * WORK, int LWORK, int * INFO) const {
  sgehrd_(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::GEHRD(int N, int ILO, int IHI, double * A, int LDA, double * TAU, 
			 double * WORK, int LWORK, int * INFO) const {
  dgehrd_(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, float * H, int LDH, 
			  float * WR, float * WI, float * Z, int LDZ, float * WORK, int LWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  shseqr_(JOB_MACRO, COMPZ_MACRO, &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::HSEQR( char JOB, char COMPZ, int N, int ILO, int IHI, double * H, int LDH, 
			  double * WR, double * WI, double * Z, int LDZ, double * WORK, int LWORK, 
			  int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dhseqr_(JOB_MACRO, COMPZ_MACRO, &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::ORGHR( int N, int ILO, int IHI, float * A, int LDA, float * TAU, 
			  float * WORK, int LWORK, int * INFO) const {
  sorghr_( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::ORGHR( int N, int ILO, int IHI, double * A, int LDA, double * TAU, 
			  double * WORK, int LWORK, int * INFO) const {
  dorghr_( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, float * A, int LDA, 
			  float * TAU, float * C, int LDC, float * WORK, int LWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  sormhr_(SIDE_MACRO, TRANS_MACRO, &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::ORMHR( char SIDE, char TRANS, int M, int N, int ILO, int IHI, double * A, int LDA, 
			  double * TAU, double * C, int LDC, double * WORK, int LWORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dormhr_(SIDE_MACRO, TRANS_MACRO, &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void Petra_LAPACK::TREVC( char SIDE, char HOWMNY, int * SELECT, int N, float * T, int LDT, float *VL, int LDVL,
			  float * VR, int LDVR, int MM, int * M, float * WORK, int * INFO) const {

#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  strevc_(SIDE_MACRO, HOWMNY_MACRO, SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void Petra_LAPACK::TREVC( char SIDE, char HOWMNY, int * SELECT, int N, double * T, int LDT, double *VL, int LDVL,
			  double * VR, int LDVR, int MM, int * M, double * WORK, int * INFO) const {

#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  dtrevc_(SIDE_MACRO, HOWMNY_MACRO, SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void Petra_LAPACK::TREXC( char COMPQ, int N, float * T, int LDT, float * Q, int LDQ, int IFST, int ILST, 
			  float * WORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  strexc_( COMPQ_MACRO, &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
void Petra_LAPACK::TREXC( char COMPQ, int N, double * T, int LDT, double * Q, int LDQ, int IFST, int ILST, 
			  double * WORK, int * INFO) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  dtrexc_( COMPQ_MACRO, &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
float Petra_LAPACK::SLAMCH( char CMACH) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  return(slamch_( CMACH_MACRO));
}
//=============================================================================
double Petra_LAPACK::DLAMCH( char CMACH) const {
#if defined (INTEL_CXML)
	unsigned int one=1;
#endif
  return(dlamch_( CMACH_MACRO));
}
  


