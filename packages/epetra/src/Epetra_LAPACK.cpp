
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
void Epetra_LAPACK::POTRF( const char UPLO, const int N, float * A, const int LDA, int * INFO) const {
  SPOTRF_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRF( const char UPLO, const int N, double * A, const int LDA, int * INFO) const {
  DPOTRF_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRS( const char UPLO, const int N, const int NRHS, const float * A, const int LDA, 
			  float * X, const int LDX, int * INFO) const {
  SPOTRS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRS( const char UPLO, const int N, const int NRHS, const double * A, const int LDA, 
			  double * X, const int LDX, int * INFO) const {
  DPOTRS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRI( const char UPLO, const int N, float * A, const int LDA, int * INFO) const {
  SPOTRI_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POTRI( const char UPLO, const int N, double * A, const int LDA, int * INFO) const {
  DPOTRI_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
}
//=============================================================================
void Epetra_LAPACK::POCON( const char UPLO, const int N, const float * A, const int LDA, const float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
  SPOCON_F77(CHAR_MACRO(UPLO), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POCON( const char UPLO, const int N, const double * A, const int LDA, const double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
  DPOCON_F77(CHAR_MACRO(UPLO), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSV( const char UPLO, const int N, const int NRHS, float * A, const int LDA, 
			  float * X, const int LDX, int * INFO) const {
  SPOSV_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSV( const char UPLO, const int N, const int NRHS, double * A, const int LDA, 
			  double * X, const int LDX, int * INFO) const {
  DPOSV_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POEQU(const int N, const float * A, const int LDA, float * S, float * SCOND, 
			 float * AMAX, int * INFO) const {
  SPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::POEQU(const int N, const double * A, const int LDA, double * S, double * SCOND,
			double * AMAX, int * INFO) const {		 
  DPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::PORFS(const char UPLO, const int N, const int NRHS, const float * A, const int LDA, const float * AF, const int LDAF, 
	     const float * B, const int LDB, float * X, const int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SPORFS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::PORFS(const char UPLO, const int N, const int NRHS, const double * A, const int LDA, const double * AF, const int LDAF, 
	     const double * B, const int LDB, double * X, const int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
  DPORFS_F77( CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF,B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSVX(const char FACT, const char UPLO, const int N, const int NRHS, float * A, const int LDA, float * AF, const int LDAF, 
	     const char EQUED, float * S, float * B, const int LDB, float * X, const int LDX, float * RCOND, 
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, CHAR_MACRO(EQUED), S, B, &LDB, X, &LDX, 
	  RCOND, FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::POSVX(const char FACT, const char UPLO, const int N, const int NRHS, double * A, const int LDA, double * AF, const int LDAF, 
	     const char EQUED, double * S, double * B, const int LDB, double * X, const int LDX, double * RCOND, 
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
  DPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, CHAR_MACRO(EQUED), S, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}

// General linear systems
//=============================================================================
void Epetra_LAPACK::GELS( const char TRANS, const int M, const int N, const int NRHS,
			  double* A, const int LDA, double* B, const int LDB, double* WORK, const int LWORK, 
			  int * INFO) const {
	DGELS_F77 (CHAR_MACRO(TRANS), &M, &N, &NRHS, A, &LDA, B, &LDB, WORK, &LWORK, INFO );
				
}
//=============================================================================
void Epetra_LAPACK::GETRF( const int M, const int N, float * A, const int LDA, int * IPIV, int * INFO) const {
  SGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRF( const int M, const int N, double * A, const int LDA, int * IPIV, int * INFO) const {
  DGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEQRF( const int M, const int N,  float * A, const int LDA,  float * TAU,  float * WORK, const int LWORK, int * INFO) const {
  SGEQRF_F77(&M, &N, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEQRF( const int M, const int N, double * A, const int LDA, double * TAU, double * WORK, const int LWORK, int * INFO) const {
  DGEQRF_F77(&M, &N, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESVD(const char JOBU, const char JOBVT, const int M, const int N, float * A, 
			  const int LDA, float * S, float * U,
			  const int LDU, float * VT, const int LDVT, float * WORK, 
			  const int * LWORK, int * INFO) const {
  SGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &M, &N, A, &LDA, S, U, &LDU,
	     VT, &LDVT, WORK, LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESVD(const char JOBU, const char JOBVT, const int M, const int N, double * A, 
			  const int LDA, double * S, double * U,
			  const int LDU, double * VT, const int LDVT, double * WORK, 
			  const int * LWORK, int * INFO) const {
  DGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &M, &N, A, &LDA, S, U, &LDU,
	     VT, &LDVT, WORK, LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRS( const char TRANS, const int N, const int NRHS, const float * A, const int LDA, 
			  const int * IPIV, float * X, const int LDX, int * INFO) const {
  SGETRS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRS( const char TRANS, const int N, const int NRHS, const double * A, const int LDA, 
			  const int * IPIV, double * X, const int LDX, int * INFO) const {
  DGETRS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRI( const int N, float * A, const int LDA, int * IPIV, 
			  float * WORK, const int * LWORK, int * INFO) const {
  SGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GETRI( const int N, double * A, const int LDA, int * IPIV, 
			  double * WORK, const int * LWORK, int * INFO) const {
  DGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GECON( const char NORM, const int N, const float  * A, const int LDA, const float ANORM, 
			  float * RCOND, float * WORK, int * IWORK, 
			  int * INFO) const {
  SGECON_F77(CHAR_MACRO(NORM), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GECON( const char NORM, const int N, const double * A, const int LDA, const double ANORM, 
			  double * RCOND, double * WORK, int * IWORK, 
			  int * INFO) const {
  DGECON_F77(CHAR_MACRO(NORM), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESV( const int N, const int NRHS, float * A, const int LDA, int * IPIV, 
			  float * X, const int LDX, int * INFO) const {
  SGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESV( const int N, const int NRHS, double * A, const int LDA, int * IPIV, 
			  double * X, const int LDX, int * INFO) const {
  DGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEEQU(const int M, const int N, const float * A, const int LDA, float * R, float * C, 
			 float * ROWCND, float * COLCND, float * AMAX, int * INFO) const {
  SGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEEQU(const int M, const int N, const double * A, const int LDA, double * R, double * C,  
			 double * ROWCND, double * COLCND, double * AMAX, int * INFO) const {
  DGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}
//=============================================================================
void Epetra_LAPACK::GERFS(const char TRANS, const int N, const int NRHS, const float * A, const int LDA, const float * AF, const int LDAF, 
	      const int * IPIV, const float * B, const int LDB, float * X, const int LDX,
	     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SGERFS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GERFS(const char TRANS, const int N, const int NRHS, const double * A, const int LDA, const double * AF, const int LDAF, 
	     const int * IPIV, const double * B, const int LDB, double * X, const int LDX,
	     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
  DGERFS_F77( CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX,
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESVX(const char FACT, const char TRANS, const int N, const int NRHS, float * A, const int LDA, float * AF, const int LDAF, 
			 int * IPIV, const char EQUED, float * R, float * C, float * B, const int LDB, float * X, const int LDX, float * RCOND, 
			 float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const {
  SGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, CHAR_MACRO(EQUED), R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESVX(const char FACT, const char TRANS, const int N, const int NRHS, double * A, const int LDA, double * AF, const int LDAF, 
			 int * IPIV, const char EQUED, double * R, double * C, double * B, const int LDB, double * X, const int LDX, double * RCOND, 
			 double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const {
 DGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, CHAR_MACRO(EQUED), R, C, B, &LDB, X, &LDX, RCOND, 
	  FERR, BERR, WORK, IWORK, INFO);
}




//=============================================================================
void Epetra_LAPACK::GEHRD(const int N, const int ILO, const int IHI, float * A, const int LDA, float * TAU, 
			 float * WORK, const int LWORK, int * INFO) const {
  SGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEHRD(const int N, const int ILO, const int IHI, double * A, const int LDA, double * TAU, 
			 double * WORK, const int LWORK, int * INFO) const {
  DGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::HSEQR( const char JOB, const char COMPZ, const int N, const int ILO, const int IHI, float * H, const int LDH, 
			  float * WR, float * WI, float * Z, const int LDZ, float * WORK, const int LWORK, 
			  int * INFO) const {
  SHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::HSEQR( const char JOB, const char COMPZ, const int N, const int ILO, const int IHI, double * H, const int LDH, 
			  double * WR, double * WI, double * Z, const int LDZ, double * WORK, const int LWORK, 
			  int * INFO) const {
  DHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORGQR( const int M, const int N, const int K, float * A, const int LDA, float * TAU, 
			  float * WORK, const int LWORK, int * INFO) const {
  SORGQR_F77( &M, &N, &K, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORGQR( const int M, const int N, const int K, double * A, const int LDA, double * TAU, 
			  double * WORK, const int LWORK, int * INFO) const {
  DORGQR_F77( &M, &N, &K, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORGHR( const int N, const int ILO, const int IHI, float * A, const int LDA, float * TAU, 
			  float * WORK, const int LWORK, int * INFO) const {
  SORGHR_F77( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORGHR( const int N, const int ILO, const int IHI, double * A, const int LDA, double * TAU, 
			  double * WORK, const int LWORK, int * INFO) const {
  DORGHR_F77( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORMHR( const char SIDE, const char TRANS, const int M, const int N, const int ILO, const int IHI, const float * A, const int LDA, 
			  const float * TAU, float * C, const int LDC, float * WORK, const int LWORK, int * INFO) const {
  SORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::ORMHR( const char SIDE, const char TRANS, const int M, const int N, const int ILO, const int IHI, const double * A, const int LDA, 
			  const double * TAU, double * C, const int LDC, double * WORK, const int LWORK, int * INFO) const {
  DORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::LARFT( const char DIRECT, const char STOREV, const int N, const int K, float * V, const int LDV, float * TAU, float * T, const int LDT) const {
  SLARFT_F77(CHAR_MACRO(DIRECT), CHAR_MACRO(STOREV), &N, &K, V, &LDV, TAU, T, &LDT);
}
//=============================================================================
void Epetra_LAPACK::LARFT( const char DIRECT, const char STOREV, const int N, const int K, double * V, const int LDV, double * TAU, double * T, const int LDT) const {
  DLARFT_F77(CHAR_MACRO(DIRECT), CHAR_MACRO(STOREV), &N, &K, V, &LDV, TAU, T, &LDT);
}
//=============================================================================
void Epetra_LAPACK::TREVC( const char SIDE, const char HOWMNY, int * SELECT, const int N, const float * T, const int LDT, float *VL, const int LDVL,
			  float * VR, const int LDVR, const int MM, int * M, float * WORK, int * INFO) const {

  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  STREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::TREVC( const char SIDE, const char HOWMNY, int * SELECT, const int N, const double * T, const int LDT, double *VL, const int LDVL,
			  double * VR, const int LDVR, const int MM, int * M, double * WORK, int * INFO) const {

  if (HOWMNY=='S') *INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)

  else  DTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::TREXC( const char COMPQ, const int N, float * T, const int LDT, float * Q, const int LDQ, int IFST, int ILST, 
			  float * WORK, int * INFO) const {
  STREXC_F77( CHAR_MACRO(COMPQ), &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::TREXC( const char COMPQ, const int N, double * T, const int LDT, double * Q, const int LDQ, int IFST, int ILST, 
			  double * WORK, int * INFO) const {
  DTREXC_F77( CHAR_MACRO(COMPQ), &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::LAMCH( const char CMACH, float & T) const {
  T = SLAMCH_F77( CHAR_MACRO(CMACH));
}
//=============================================================================
void Epetra_LAPACK::LAMCH( const char CMACH, double & T) const {
  T = DLAMCH_F77( CHAR_MACRO(CMACH));
}
//=============================================================================
void Epetra_LAPACK::GGSVD(const char JOBU, const char JOBV, const char JOBQ, const int M, const int N, const int P, int * K, int * L,  
			  double* A,  const int LDA,  double* B,  const int LDB,
                          double* ALPHA,  double* BETA,  double* U,  const int LDU, double* V, const int LDV, double* Q, const int LDQ, double* WORK, int* IWORK,
                          int* INFO) const {
  DGGSVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBV), CHAR_MACRO(JOBQ), &M, &N, &P, K, L,  A,  &LDA,  B,  &LDB,
	     ALPHA,  BETA,  U,  &LDU, V, &LDV, Q, &LDQ, WORK, IWORK, INFO);
}
  
//=============================================================================
void Epetra_LAPACK::GGSVD(const char JOBU, const char JOBV, const char JOBQ, const int M, const int N, const int P, int * K, int * L,  
			  float* A,  const int LDA,  float* B,  const int LDB,
                          float* ALPHA,  float* BETA,  float* U,  const int LDU, float* V, const int LDV, float* Q, const int LDQ, float* WORK, int* IWORK,
                          int* INFO) const {
  SGGSVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBV), CHAR_MACRO(JOBQ), &M, &N, &P, K, L,  A,  &LDA,  B,  &LDB,
	     ALPHA,  BETA,  U,  &LDU, V, &LDV, Q, &LDQ, WORK, IWORK, INFO);
}
  
//=============================================================================
void Epetra_LAPACK::GEEV(const char JOBVL, const char JOBVR, const int N, double* A, const int LDA, double* WR, double* WI, 
			 double* VL, const int LDVL, double* VR, const int LDVR, double* WORK, const int LWORK, int* INFO) const {

  DGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &N, A, &LDA, WR, WI, VL, &LDVL, VR,  &LDVR,
	    WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEEV(const char JOBVL, const char JOBVR, const int N, float* A, const int LDA, float* WR, float* WI, 
			 float* VL, const int LDVL, float* VR, const int LDVR, float* WORK, const int LWORK, int* INFO) const {

  SGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &N, A, &LDA, WR, WI, VL, &LDVL, VR,  &LDVR,
	    WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SPEV(const char JOBZ, const char UPLO, const int N, double* AP, double* W, double* Z, int LDZ, double* WORK, int* INFO) const {

  DSPEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, AP, W, Z, &LDZ, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SPEV(const char JOBZ, const char UPLO, const int N, float* AP, float* W, float* Z, int LDZ, float* WORK, int* INFO) const {

  SSPEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, AP, W, Z, &LDZ, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SPGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, double* AP, double* BP, double* W, double* Z, const int LDZ, double* WORK, int* INFO) const {

  DSPGV_F77(&ITYPE, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, AP, BP, W, Z, &LDZ, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SPGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, float* AP, float* BP, float* W, float* Z, const int LDZ, float* WORK, int* INFO) const {

  SSPGV_F77(&ITYPE, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, AP, BP, W, Z, &LDZ, WORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SYEV(const char JOBZ, const char UPLO, const int N, double* A, const int LDA, double* W, double* WORK, const int LWORK, int* INFO) const{

  DSYEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &LDA, W, WORK, &LWORK, INFO);

}
//=============================================================================
void Epetra_LAPACK::SYEV(const char JOBZ, const char UPLO, const int N, float* A, const int LDA, float* W, float* WORK, const int LWORK, int* INFO) const{

  SSYEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &LDA, W, WORK, &LWORK, INFO);

}
//=============================================================================
void Epetra_LAPACK::SYEVD(const char JOBZ, const char UPLO,  const int N,  double* A,  const int LDA,  double* W,  double* WORK,  
			  const int LWORK,  int* IWORK, const int LIWORK, int* INFO) const {

  DSYEVD_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SYEVD(const char JOBZ, const char UPLO,  const int N,  float* A,  const int LDA,  float* W,  float* WORK,  
			  const int LWORK,  int* IWORK, const int LIWORK, int* INFO) const {

  SSYEVD_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SYEVX(const char JOBZ, const char RANGE, const char UPLO,  const int N,  double* A,  const int LDA,  
			   const double* VL,  const double* VU,  const int* IL,  const int* IU,
			   const double ABSTOL,  int * M,  double* W,  double* Z,  const int LDZ, double* WORK, 
			  const int LWORK, int* IWORK, int* IFAIL,
			  int* INFO) const {

  DSYEVX_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(RANGE), CHAR_MACRO(UPLO), &N,  A,  &LDA,  VL,  VU,  IL,  IU,
	 &ABSTOL,  M,  W,  Z,  &LDZ, WORK, &LWORK, IWORK, IFAIL, INFO);
}
//=============================================================================
void Epetra_LAPACK::SYEVX(const char JOBZ, const char RANGE, const char UPLO,  const int N,  float* A,  const int LDA,  
			   const float* VL,  const float* VU,  const int* IL,  const int* IU,
			   const float ABSTOL,  int * M,  float* W,  float* Z,  const int LDZ, float* WORK, 
			  const int LWORK, int* IWORK, int* IFAIL,
			  int* INFO) const {

  SSYEVX_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(RANGE), CHAR_MACRO(UPLO), &N,  A,  &LDA,  VL,  VU,  IL,  IU,
	 &ABSTOL,  M,  W,  Z,  &LDZ, WORK, &LWORK, IWORK, IFAIL, INFO);
}
//=============================================================================
void Epetra_LAPACK::SYGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, double* A, 
			 const int LDA, double* B, const int LDB, double* W, double* WORK, const int LWORK, int* INFO) const{

  DSYGV_F77(&ITYPE, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &LDA, B, &LDB, W, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SYGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, float* A, 
			 const int LDA, float* B, const int LDB, float* W, float* WORK, const int LWORK, int* INFO) const{

  SSYGV_F77(&ITYPE, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &LDA, B, &LDB, W, WORK, &LWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::SYGVX(const int ITYPE, const char JOBZ, const char RANGE, const char UPLO, const int N, 
	     double* A, const int LDA, double* B, const int LDB, const double* VL, const double* VU,
	     const int* IL, const int* IU, const double ABSTOL, int* M, double* W, double* Z, 
	     const int LDZ,  double* WORK,  const int LWORK,  int* IWORK,
			  int* IFAIL, int* INFO) const {

#ifdef EPETRA_LAPACK3
  DSYGVX_F77(&ITYPE, CHAR_MACRO(JOBZ), CHAR_MACRO(RANGE), CHAR_MACRO(UPLO), &N,  A,  &LDA, B, &LDB, VL, VU,
	     IL, IU, &ABSTOL, M, W, Z, &LDZ,  WORK,  &LWORK,  IWORK,
	     IFAIL, INFO);
#else

  Epetra_Object obj;
  obj.ReportError("SYGVX requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::SYGVX(const int ITYPE, const char JOBZ, const char RANGE, const char UPLO, const int N, 
	     float* A, const int LDA, float* B, const int LDB, const float* VL, const float* VU,
	     const int* IL, const int* IU, const float ABSTOL, int* M, float* W, float* Z, 
	     const int LDZ,  float* WORK,  const int LWORK,  int* IWORK,
			  int* IFAIL, int* INFO) const {

#ifdef EPETRA_LAPACK3
  SSYGVX_F77(&ITYPE, CHAR_MACRO(JOBZ), CHAR_MACRO(RANGE), CHAR_MACRO(UPLO), &N,  A,  &LDA, B, &LDB, VL, VU,
                          IL, IU, &ABSTOL, M, W, Z, &LDZ,  WORK,  &LWORK,  IWORK,
	     IFAIL, INFO);
#else

  Epetra_Object obj;
  obj.ReportError("SYGVX requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::SYEVR(const char JOBZ, const char RANGE, const char UPLO,  const int N,  double* A,  const int LDA,  
			  const double* VL,  const double* VU,  const int *IL,  const int *IU,
                          const double ABSTOL,  int* M,  double* W,  double* Z, const int LDZ, int* ISUPPZ, double* WORK, const int LWORK, int* IWORK,
                          const int LIWORK, int* INFO) const { 

#ifdef EPETRA_LAPACK3
  DSYEVR_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(RANGE), CHAR_MACRO(UPLO), &N,  A,  &LDA, VL,  VU,  IL,  IU,
	     &ABSTOL,  M,  W,  Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK,
	     &LIWORK, INFO);
#else
  Epetra_Object obj;
  obj.ReportError("SYEVR requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::SYEVR(const char JOBZ, const char RANGE, const char UPLO,  const int N,  float* A,  const int LDA,  
			  const float* VL,  const float* VU,  const int *IL,  const int *IU,
                          const float ABSTOL,  int* M,  float* W,  float* Z, const int LDZ, int* ISUPPZ, float* WORK, const int LWORK, int* IWORK,
                          const int LIWORK, int* INFO) const { 

#ifdef EPETRA_LAPACK3
  SSYEVR_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(RANGE), CHAR_MACRO(UPLO), &N,  A,  &LDA, VL,  VU,  IL,  IU,
	     &ABSTOL,  M,  W,  Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK,
	     &LIWORK, INFO);
#else
  Epetra_Object obj;
  obj.ReportError("SYEVR requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int N, double* A, const int LDA, double* WR, double* WI,  double* VL,
	     const int LDVL,  double* VR,  const int LDVR,  int* ILO,  int* IHI,  double* SCALE, double* ABNRM, double* RCONDE,
			  double* RCONDV, double* WORK, const int LWORK, int* IWORK, int* INFO) const{

  DGEEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &N, A, &LDA, WR, WI,  VL,
	     &LDVL,  VR,  &LDVR,  ILO,  IHI,  SCALE, ABNRM, RCONDE,
	     RCONDV, WORK, &LWORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int N, float* A, const int LDA, float* WR, float* WI,  float* VL,
	     const int LDVL,  float* VR,  const int LDVR,  int* ILO,  int* IHI,  float* SCALE, float* ABNRM, float* RCONDE,
			  float* RCONDV, float* WORK, const int LWORK, int* IWORK, int* INFO) const{

  SGEEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &N, A, &LDA, WR, WI,  VL,
	     &LDVL,  VR,  &LDVR,  ILO,  IHI,  SCALE, ABNRM, RCONDE,
	     RCONDV, WORK, &LWORK, IWORK, INFO);
}
//=============================================================================
void Epetra_LAPACK::GESDD(const char JOBZ, const int M, const int N, double* A, const int LDA,  double* S,  
			  double* U,  const int LDU,  double* VT,  const int LDVT,  double* WORK,
			  const int LWORK, int* IWORK, int* INFO) const{

#ifdef EPETRA_LAPACK3
  DGESDD_F77(CHAR_MACRO(JOBZ), &M, &N, A, &LDA,  S,  U,  &LDU,  VT,  &LDVT,  WORK,
	     &LWORK, IWORK, INFO);
#else
  Epetra_Object obj;
  obj.ReportError("GESDD requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::GESDD(const char JOBZ, const int M, const int N, float* A, const int LDA,  float* S,  
			  float* U,  const int LDU,  float* VT,  const int LDVT,  float* WORK,
			  const int LWORK, int* IWORK, int* INFO) const{

#ifdef EPETRA_LAPACK3
  SGESDD_F77(CHAR_MACRO(JOBZ), &M, &N, A, &LDA,  S,  U,  &LDU,  VT,  &LDVT,  WORK,
	     &LWORK, IWORK, INFO);
#else
  Epetra_Object obj;
  obj.ReportError("GESDD requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::GGEV(const char JOBVL,  const char JOBVR,  const int N,  double* A,  
			 const int LDA,  double* B, const int LDB, double* ALPHAR, double* ALPHAI,
			 double* BETA, double* VL, const int LDVL, double* VR, const int 
			 LDVR, double* WORK, const int LWORK, int* INFO) const{

#ifdef EPETRA_LAPACK3
  DGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &N,  A,  &LDA,  B, &LDB, ALPHAR, ALPHAI,
	    BETA, VL, &LDVL, VR, &LDVR, WORK, &LWORK, INFO);
#else
  Epetra_Object obj;
  obj.ReportError("GGEV requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::GGEV(const char JOBVL,  const char JOBVR,  const int N,  float* A,  
			 const int LDA,  float* B, const int LDB, float* ALPHAR, float* ALPHAI,
			 float* BETA, float* VL, const int LDVL, float* VR, const int 
			 LDVR, float* WORK, const int LWORK, int* INFO) const {

#ifdef EPETRA_LAPACK3
  SGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &N,  A,  &LDA,  B, &LDB, ALPHAR, ALPHAI,
	    BETA, VL, &LDVL, VR, &LDVR, WORK, &LWORK, INFO);
#else
  Epetra_Object obj;
  obj.ReportError("GGEV requires LAPACK Version 3.  Compile Epetra with -DEPETRA_LAPACK3 and link with LAPACK 3 library", -1);
#endif
}
//=============================================================================
void Epetra_LAPACK::GGLSE(const int M, const int N, const int P, double* A, const int LDA, double* B, const int LDB, 
			  double* C, double* D, double* X, double* WORK, const int LWORK, int* INFO) const{
  DGGLSE_F77(&M, &N, &P, A, &LDA, B, &LDB, C, D, X, WORK, &LWORK,  INFO);
}
//=============================================================================
void Epetra_LAPACK::GGLSE(const int M, const int N, const int P, float* A, const int LDA, float* B, const int LDB, 
			  float* C, float* D, float* X, float* WORK, const int LWORK, int* INFO) const{
  SGGLSE_F77(&M, &N, &P, A, &LDA, B, &LDB, C, D, X, WORK, &LWORK,  INFO);
}
//=============================================================================
void Epetra_LAPACK::TRTRS(const char UPLO, const char TRANS, const char DIAG, const int N, const int NRHS, const float *A,
                          const int LDA, float *B, const int LDB, int *INFO) const{
  STRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &N, &NRHS, A, &LDA, B, &LDB, INFO);
}
//=============================================================================
void Epetra_LAPACK::TRTRS(const char UPLO, const char TRANS, const char DIAG, const int N, const int NRHS, const double *A,
                          const int LDA, double *B, const int LDB, int *INFO) const{
  DTRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &N, &NRHS, A, &LDA, B, &LDB, INFO);
}

