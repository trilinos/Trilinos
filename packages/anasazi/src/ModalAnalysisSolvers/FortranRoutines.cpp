//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the report
//
// " A comparison of algorithms for modal analysis in the absence 
//   of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//  Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://software.sandia.gov/trilinos/ ).
//
// The distribution of this software follows also the rules defined in Trilinos.
// This notice shall be marked on any reproduction of this software, in whole or
// in part.
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

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

#include "FortranRoutines.h"

// Double precision BLAS 1 //

void FortranRoutines::SCAL_INCX(int N, double ALPHA, double *X, int incX) const {
  DSCAL_F77(&N, &ALPHA, X, &incX);
  return;
}

void FortranRoutines::SWAP(int N, double *X, int incx, double *Y, int incy) const {
  F77_FUNC(dswap,DSWAP)(&N, X, &incx, Y, &incy);
  return;
}


// Double precision LAPACK //


void FortranRoutines::GEQRF(int M, int N, double *A, int lda, double *tau, double *work,
                            int lwork, int *info) const {
  F77_FUNC(dgeqrf,DGEQRF)(&M, &N, A, &lda, tau, work, &lwork, info);
  return;
}

void FortranRoutines::ORMQR(char SIDE, char TRANS, int M, int N, int K, double *A, int lda, 
                            double *tau, double *C, int ldc, double *work, int lwork,
                            int *info) const {
  F77_FUNC(dormqr,DORMQR)(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &K, A, &lda, tau, 
                          C, &ldc, work, &lwork, info);
  return;
}

void FortranRoutines::SPEV(char JOBZ, char UPLO, int N, double *A, double *W, double *Z,
                           int ldz, double *work, int *info) const {
  F77_FUNC(dspev,DSPEV)(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, W, Z, &ldz, work, info);
  return;
}

void FortranRoutines::STEQR(char COMPZ, int N, double *D, double *E, double *Z, int ldz, 
                            double *work, int *info) const {
  F77_FUNC(dsteqr,DSTEQR)(CHAR_MACRO(COMPZ), &N, D, E, Z, &ldz, work, info);
  return;
}

void FortranRoutines::SYEV(char JOBZ, char UPLO, int N, double *A, int lda, double *W,
                           double *work, int lwork, int *info) const {
  F77_FUNC(dsyev,DSYEV)(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &lda, W, work, &lwork, info);
  return;
}

void FortranRoutines::SYGV(int itype, char JOBZ, char UPLO, int N, double *A, int lda, 
                           double *B, int ldb, double *W, double *work, int lwork, 
                           int *info) const {
  F77_FUNC(dsygv,DSYGV)(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &N, A, &lda, B, &ldb,
           W, work, &lwork, info);
  return;
}

int FortranRoutines::LAENV(int ispec, char *NAME, char *OPTS, int N1, int N2, int N3,
                           int N4, int len_name, int len_opts) const {
#if defined (INTEL_CXML)
  return F77_FUNC(ilaenv,ILAENV)(&ispec, NAME, len_name, OPTS, len_opts, &N1, &N2, &N3, &N4);
#else
  return F77_FUNC(ilaenv,ILAENV)(&ispec, NAME, OPTS, &N1, &N2, &N3, &N4, len_name, len_opts);
#endif
}


// Double precision ARPACK routines


void FortranRoutines::SAUPD(int *ido, char BMAT, int N, char *which, int nev, double tol, 
                            double *resid, int ncv, double *V, int ldv, int *iparam, 
                            int *ipntr, double *workd, double *workl, int lworkl, int *info,
                            int verbose) const {
#if defined (INTEL_CXML)
  F77_FUNC(mydsaupd,MYDSAUPD)(ido, &BMAT, 1, &N, which, 2, &nev, &tol, resid, &ncv, V, &ldv,
           iparam, ipntr, workd, workl, &lworkl, info, &verbose);
#else
  F77_FUNC(mydsaupd,MYDSAUPD)(ido, &BMAT, &N, which, &nev, &tol, resid, &ncv, V, &ldv,
           iparam, ipntr, workd, workl, &lworkl, info, &verbose, 1, 2);
#endif
  return;
}

void FortranRoutines::SEUPD(LOGICAL rvec, char HOWMNY, LOGICAL *select, double *D, 
                            double *Z, int ldz, double sigma, char BMAT, int N, 
                            char *which, int nev, double tol, double *resid, int ncv, double *V,
                            int ldv, int *iparam, int *ipntr, double *workd, double *workl,
                            int lworkl, int *info) const {
#if defined (INTEL_CXML)
  F77_FUNC(dseupd,DSEUPD)(&rvec, &HOWMNY, 1, select, D, Z, &ldz, &sigma, &BMAT, 1, &N,
           which, 2, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl,
           info);
#else
  F77_FUNC(dseupd,DSEUPD)(&rvec, &HOWMNY, select, D, Z, &ldz, &sigma, &BMAT, &N, which, 
           &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, info, 
           1, 1, 2);
#endif
  return;
}

#ifdef EPETRA_MPI

// Double precision PARPACK routines

void FortranRoutines::PSAUPD(MPI_Comm MyComm, int *ido, char BMAT, int N, char *which, int nev, 
                             double tol, double *resid, int ncv, double *V, int ldv, int *iparam,
                             int *ipntr, double *workd, double *workl, int lworkl, int *info, 
                             int verbose) const {
#if defined (INTEL_CXML)
  F77_FUNC(mypdsaupd,MYPDSAUPD)(&MyComm, ido, &BMAT, 1, &N, which, 2, &nev, &tol, resid, &ncv, 
           V, &ldv, iparam, ipntr, workd, workl, &lworkl, info, &verbose);
#else
  F77_FUNC(mypdsaupd,MYPDSAUPD)(&MyComm, ido, &BMAT, &N, which, &nev, &tol, resid, &ncv, V, &ldv,
           iparam, ipntr, workd, workl, &lworkl, info, &verbose, 1, 2);
#endif
  return;
}

void FortranRoutines::PSEUPD(MPI_Comm MyComm, LOGICAL rvec, char HOWMNY, LOGICAL *select, 
                             double *D, double *Z, int ldz, double sigma, char BMAT, int N,
                             char *which, int nev, double tol, double *resid, int ncv, double *V,
                             int ldv, int *iparam, int *ipntr, double *workd, double *workl,
                             int lworkl, int *info) const {
#if defined (INTEL_CXML)
  F77_FUNC(pdseupd,PDSEUPD)(&MyComm, &rvec, &HOWMNY, 1, select, D, Z, &ldz, &sigma, &BMAT, 1, &N,
           which, 2, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, info);
#else
  F77_FUNC(pdseupd,PDSEUPD)(&MyComm, &rvec, &HOWMNY, select, D, Z, &ldz, &sigma, &BMAT, &N,
           which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, info, 
           1, 1, 2);
#endif
  return;
}

#endif

