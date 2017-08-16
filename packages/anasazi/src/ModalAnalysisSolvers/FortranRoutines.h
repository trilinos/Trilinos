// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ***********************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#ifndef FORTRAN_ROUTINES
#define FORTRAN_ROUTINES

#include "Epetra_BLAS_wrappers.h"
#include "Epetra_LAPACK_wrappers.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#endif

// LOGICAL as 4 bytes
typedef int LOGICAL;

#ifdef __cplusplus
extern "C" {
#endif

// Double precision BLAS 1 //
void PREFIX F77_FUNC(dswap,DSWAP)(int *n, double x[], int* incx, double y[], int* incy);

// Double precision LAPACK //
void PREFIX F77_FUNC(dgeqrf,DGEQRF)(int *M, int *N, double *A, int *lda, double *tau, 
                     double *work, int *lwork, int *info);
void PREFIX F77_FUNC(dormqr,DORMQR)(Epetra_fcd, Epetra_fcd, int *M, int *N, int *K, double *A,
                     int *lda, double *tau, double *C, int *ldc, double *work, int *lwork,
                     int *info);
void PREFIX F77_FUNC(dsteqr,DSTEQR)(Epetra_fcd, int *N, double *D, double *E, double *Z,
                     int *ldz, double *work, int *info);

#if defined (INTEL_CXML)
int PREFIX F77_FUNC(ilaenv,ILAENV)(int *ispec, char *NAME, unsigned int len_name, char *OPTS, 
                    unsigned int len_opts, int *N1, int *N2, int *N3, int *N4);
#else
int PREFIX F77_FUNC(ilaenv,ILAENV)(int *ispec, char *NAME, char *OPTS, int *N1, int *N2,
                    int *N3, int *N4, int len_name, int len_opts);
#endif

// Double precision customized ARPACK routines //
#if defined (INTEL_CXML)
void PREFIX F77_FUNC(mydsaupd,MYDSAUPD)(int *, char *, unsigned int, int *, char *, 
                     unsigned int, int *, double *, double *, int *, double *, int *, int *,
                     int *, double *, double *, int *, int *, int *);
void PREFIX F77_FUNC(dseupd,DSEUPD)(LOGICAL *rvec, char *HOWMNY, unsigned int len_howny,
                     LOGICAL *select, double *D, double *Z, int *ldz, double *sigma, char *BMAT,
                     unsigned int len_bmat, int *N, char *which, unsigned int len_which, 
                     int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv,
                     int *iparam, int *ipntr, double *workd, double *workl, int *lworkl,
                     int *info);
#else
void PREFIX F77_FUNC(mydsaupd,MYDSAUPD)(int *, char *, int *, char *, int *, 
                     double *, double *, int *, double *, int *, int *, int *, double *,
                     double *, int *, int *, int *, int, int);
void PREFIX F77_FUNC(dseupd,DSEUPD)(LOGICAL *rvec, char *HOWMNY, LOGICAL *select, double *D, 
                     double *Z, int *ldz, double *sigma, char *BMAT, int *N, char *which,
                     int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv,
                     int *iparam, int *ipntr, double *workd, double *workl, int *lworkl,
                     int *info, int len_howmny, int len_bmat, int len_which);
#endif

#ifdef EPETRA_MPI

// Double precision customized PARPACK routines //
#if defined (INTEL_CXML)
void PREFIX F77_FUNC(mypdsaupd,MYPDSAUPD)(MPI_Comm *, int *, char *, unsigned int, int *, 
                     char *, unsigned int, int *, double *, double *, int *, double *,
                     int *, int *, int *, double *, double *, int *, int *, int *);
void PREFIX F77_FUNC(pdseupd,PDSEUPD)(MPI_Comm *MyComm, LOGICAL *rvec, char *HOWMNY, 
                     unsigned int len_howmny, LOGICAL *select, double *D, double *Z,
                     int *ldz, double *sigma, char *BMAT, unsigned int len_bmat, int *N,
                     char *which, unsigned int len_which, int *nev, double *tol, double *resid,
                     int *ncv, double *V, int *ldv, int *iparam, int *ipntr, double *workd,
                     double *workl, int *lworkl, int *info);
#else
void PREFIX F77_FUNC(mypdsaupd,MYPDSAUPD)(MPI_Comm *, int *, char *, int *, char *, int *, 
                     double *, double *, int *, double *, int *, int *, int *, double *,
                     double *, int *, int *, int *, int, int);
void PREFIX F77_FUNC(pdseupd,PDSEUPD)(MPI_Comm *MyComm, LOGICAL *rvec, char *HOWMNY, 
                     LOGICAL *select, double *D, double *Z, int *ldz, double *sigma,
                     char *BMAT, int *N, char *which, int *nev, double *tol, double *resid, 
                     int *ncv, double *V, int *ldv, int *iparam, int *ipntr, double *workd,
                     double *workl, int *lworkl, int *info, int len_howmny,
                     int len_bmat, int len_which);
#endif

#endif

#ifdef __cplusplus
}
#endif

class FortranRoutines {

  public: 

  // Double precision BLAS 1 //
  void SCAL_INCX(int N, double ALPHA, double *X, int incX) const;
  void SWAP(int N, double *X, int incx, double *Y, int incy) const;

  // Double precision LAPACK //
  void GEQRF(int M, int N, double *A, int lda, double *tau, double *work, int lwork, 
             int *info) const;
  void ORMQR(char SIDE, char TRANS, int M, int N, int K, double *A, int lda, double *tau, 
             double *C, int ldc, double *work, int lwork, int *info) const;
  void SPEV(char JOBZ, char UPLO, int N, double *A, double *W, double *Z, int ldz, 
            double *work, int *info) const;
  void STEQR(char COMPZ, int N, double *D, double *E, double *Z, int ldz, double *work,
             int *info) const;
  void SYEV(char JOBZ, char UPLO, int N, double *A, int lda, double *W, double *work,
            int lwork, int *info) const;
  void SYGV(int itype, char JOBZ, char UPLO, int N, double *A, int lda, double *B, int ldb,
            double *W, double *work, int lwork, int *info) const;

  int LAENV(int ispec, char *NAME, char *OPTS, int N1, int N2, int N3, int N4, 
            int len_name, int len_opts) const;

  // Double precision ARPACK routines
  void SAUPD(int *ido, char BMAT, int N, char *which, int nev, double tol, double *resid,
             int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd, double *workl,
             int lworkl, int *info, int verbose) const;
  void SEUPD(LOGICAL rvec, char HOWMNY, LOGICAL *select, double *D, double *Z, int ldz,
             double sigma, char BMAT, int N, char *which, int nev, double tol, double *resid,
             int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd,
             double *workl, int lworkl, int *info) const;

#ifdef EPETRA_MPI
  // Double precision PARPACK routines
  void PSAUPD(MPI_Comm MyComm, int *ido, char BMAT, int N, char *which, int nev, double tol, 
              double *resid, int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd,
              double *workl, int lworkl, int *info, int verbose) const;
  void PSEUPD(MPI_Comm MyComm, LOGICAL rvec, char HOWMNY, LOGICAL *select, double *D, double *Z, 
              int ldz, double sigma, char BMAT, int N, char *which, int nev, double tol,
              double *resid, int ncv, double *V, int ldv, int *iparam, int *ipntr, double *workd,
              double *workl, int lworkl, int *info) const;

#endif

};

#endif
