
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

#include <stdio.h>
#include <string.h>


#if defined(CRAY_T3X)

#include "fortran.h"
#define Epetra_fcd fcd
#define PREFIX

#define dgetrf_     SGETRF
#define dgetrs_     SGETRs
#define dpotrf_     SPOTRF
#define dpotrs_     SPOTRS
#define dpotri_     SPOTRI
#define dpocon_     SPOCON
#define dposv_      SPOSV
#define dpoequ_     SPOEQU
#define dporfs_     SPORFS
#define dposvx_     SPOSVX
#define sgetrf_     SGETRF
#define sgetrs_     SGETRS
#define spotrf_     SPOTRF
#define spotrs_     SPOTRS
#define spotri_     SPOTRI
#define spocon_     SPOCON
#define sposv_      SPOSV
#define spoequ_     SPOEQU
#define sporfs_     SPORFS
#define sposvx_     SPOSVX 
#define dgels_      SGELS
#define dgeev_      SGEEV  
#define dgehrd_     SGEHRS
#define dhseqr_     SHSEQR
#define dorghr_     SORGHR
#define dormhr_     SORMHR 
#define dtrevc_     STREVC
#define dtrexc_     STREXC
#define sgels_      SGELS  
#define sgeev_      SGEEV
#define sgehrd_     SGEHRS
#define shseqr_     SHSEQR
#define sorghr_     SORGHR
#define sormhr_     SORMHR
#define strevc_     STREVC
#define strexc_     STREXC
#define slamch_     SLAMCH
#define dlamch_     DLAMCH

#elif defined(INTEL_CXML)

#define Epetra_fcd char *, unsigned int
#define PREFIX __stdcall

#define dgecon_ DGECON
#define dgeev_  DGEEV 
#define dgeequ_ DGEEQU 
#define dgehrd_ DGEHRD
#define dgels_  DGELS
#define dgerfs_ DGERFS
#define dgesv_  DGESV
#define dgesvx_ DGESVX
#define dgetrf_ DGETRF
#define dgetri_ DGETRI
#define dgetrs_ DGETRS
#define dhseqr_ DHSEQR
#define dorghr_ DORGHR
#define dormhr_ DORMHR 
#define dpotrf_ DPOTRF
#define dpotrs_ DPOTRS
#define dpotri_ DPOTRI
#define dpocon_ DPOCON
#define dposv_  DPOSV
#define dpoequ_ DPOEQU
#define dporfs_ DPORFS
#define dposvx_ DPOSVX
#define dtrevc_ DTREVC
#define dtrexc_ DTREXC

#define sgecon_ SGECON
#define sgeev_  SGEEV
#define sgeequ_ SGEEQU
#define sgehrd_ SGEHRD
#define sgels_  SGELS 
#define sgerfs_ SGERFS
#define sgesv_  SGESV 
#define sgesvx_ SGESVX
#define sgetrf_ SGETRF
#define sgetri_ SGETRI
#define sgetrs_ SGETRS
#define shseqr_ SHSEQR
#define sorghr_ SORGHR
#define sormhr_ SORMHR
#define spotrf_ SPOTRF
#define spotrs_ SPOTRS
#define spotri_ SPOTRI
#define spocon_ SPOCON
#define sposv_  SPOSV
#define spoequ_ SPOEQU
#define sporfs_ SPORFS
#define sposvx_ SPOSVX 
#define strevc_ STREVC
#define strexc_ STREXC
#define slamch_ SLAMCH
#define dlamch_ DLAMCH

#elif defined(INTEL_MKL)

#define Epetra_fcd char *
#define PREFIX

#define dgetrf_      DGETRF
#define dgetrs_      DGETRS
#define dpotrf_      DPOTRF
#define dpotrs_      DPOTRS
#define dpotri_      DPOTRI
#define dpocon_      DPOCON
#define dposv_       DPOSV
#define dpoequ_      DPOEQU
#define dporfs_      DPORFS
#define dposvx_      DPOSVX
#define sgetrf_      SGETRF
#define sgetrs_      SGETRS
#define spotrf_      SPOTRF
#define spotrs_      SPOTRS
#define spotri_      SPOTRI
#define spocon_      SPOCON
#define sposv_       SPOSV
#define spoequ_      SPOEQU
#define sporfs_      SPORFS
#define sposvx_      SPOSVX 
#define dgels_       DGELS
#define dgeev_       DGEEV  
#define dgehrd_      DGEHRS
#define dhseqr_      DHSEQR
#define dorghr_      DORGHR
#define dormhr_      DORMHR 
#define dtrevc_      DTREVC
#define dtrexc_      DTREXC
#define sgels_       SGELS  
#define sgeev_       SGEEV
#define sgehrd_      SGEHRS
#define shseqr_      SHSEQR
#define sorghr_      SORGHR
#define sormhr_      SORMHR
#define strevc_      STREVC
#define strexc_      STREXC
#define slamch_      SLAMCH
#define dlamch_      DLAMCH

#else

#define Epetra_fcd char *
#define PREFIX

#if defined(RS6000)
#define dgetrf_      dgetrf
#define dgetrs_     dgetrs
#define dpotrf_     dpotrf
#define dpotrs_     dpotrs
#define dpotri_     dpotri
#define dpocon_     dpocon
#define dposv_      dposv
#define dpoequ_     dpoequ
#define dporfs_     dporfs
#define dposvx_     dposvx
#define sgetrf_     sgetrf
#define sgetrs_     sgetrs
#define spotrf_     spotrf
#define spotrs_     spotrs
#define spotri_     spotri
#define spocon_     spocon
#define sposv_      sposv
#define spoequ_     spoequ
#define sporfs_     sporfs
#define sposvx_     sposvx
#define dgels_      dgels
#define dgeev_      dgeev
#define dgehrd_     dgehrd
#define dhseqr_     dhseqr
#define dorghr_     dorghr
#define dormhr_     dormhr
#define dtrevc_     dtrevc
#define dtrexc_     dtrexc
#define sgels_      sgels
#define sgeev_      sgeev
#define sgehrd_     sgehrd
#define shseqr_     shseqr
#define sorghr_     sorghr
#define sormhr_     sormhr
#define strevc_     strevc
#define strexc_     strexc
#define slamch_     slamch
#define dlamch_     dlamch

#endif

#endif

// Double precision LAPACK linear solvers
extern "C" void PREFIX dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info); 
extern "C" void PREFIX dgetrs_(Epetra_fcd, int* n, int* nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
extern "C" void PREFIX dgetri_(int* n, double* a, int* lda, int*ipiv, double * work , int* lwork, int* info);
extern "C" void PREFIX dgecon_(Epetra_fcd norm, int* n, double* a, int* lda, 
                       double *anorm, double * rcond, double * work,
                       int * iwork, int* info); 
extern "C" void PREFIX dgesv_(int * n, int * nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
extern "C" void PREFIX dgeequ_(int* m, int* n, double* a, int* lda, double * r, double * c, 
			double * rowcnd, double * colcnd,
                       double * amax, int* info); 
extern "C" void PREFIX dgerfs_(Epetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, int*ipiv, 
                       double * b, int * ldb, double * x, int * ldx, 
                       double * ferr, double * berr, double * work, 
                       int * iwork, int * info);
extern "C" void PREFIX dgesvx_(Epetra_fcd, Epetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, int*ipiv, Epetra_fcd, 
                       double * r, double *c, double * b, int * ldb, double * x, int * ldx, 
                       double * rcond, double * ferr, double * berr, double * work, 
                       int * iwork, int * info);

extern "C" void PREFIX dpotrf_(Epetra_fcd, int* n, double* a, int* lda, int* info); 
extern "C" void PREFIX dpotrs_(Epetra_fcd, int * n, int * nrhs, double* a,
                       int* lda, double*x , int* ldx, int* info);
extern "C" void PREFIX dpotri_(Epetra_fcd, int* n, double* a, int* lda, int* info); 
extern "C" void PREFIX dpocon_(Epetra_fcd, int* n, double* a, int* lda, 
                       double * anorm, double * rcond, double * work,
                       int * iwork, int* info); 
extern "C" void PREFIX dposv_(Epetra_fcd, int * n, int * nrhs, double* a,
                       int* lda, double*x , int* ldx, int* info);
extern "C" void PREFIX dpoequ_(int* n, double* a, int* lda, double * s, double * scond,
                       double * amax, int* info); 

extern "C" void PREFIX dporfs_(Epetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, 
                       double * b, int * ldb, double * x, int * ldx, 
                       double * ferr, double * berr, double * work, 
                       int * iwork, int * info);

extern "C" void PREFIX dposvx_(Epetra_fcd, Epetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, Epetra_fcd, 
                       double * s, double * b, int * ldb, double * x, int * ldx, 
                       double * rcond, double * ferr, double * berr, double * work, 
                       int * iwork, int * info);

// Single precision LAPACK linear solvers
extern "C" void PREFIX sgetrf_(int* m, int* n, float* a, int* lda, int* ipiv, int* info); 
extern "C" void PREFIX sgetrs_(Epetra_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
extern "C" void PREFIX sgetri_(int* n, float* a, int* lda, int*ipiv, float * work , int* lwork, int* info);
extern "C" void PREFIX sgecon_(Epetra_fcd norm, int* n, float* a, int* lda, 
                       float * anorm, float * rcond, float * work,
                       int * iwork, int* info); 
extern "C" void PREFIX sgesv_(int * n, int * nrhs, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
extern "C" void PREFIX sgeequ_(int* m, int* n, float* a, int* lda, float * r, float * c, 
			float * rowcnd, float * colcnd,
			float * amax, int* info); 
extern "C" void PREFIX sgerfs_(Epetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, int*ipiv, 
                       float * b, int * ldb, float * x, int * ldx, 
                       float * ferr, float * berr, float * work, 
                       int * iwork, int * info);
extern "C" void PREFIX sgesvx_(Epetra_fcd, Epetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, int*ipiv, Epetra_fcd, 
                       float * r, float *c, float * b, int * ldb, float * x, int * ldx, 
                       float * rcond, float * ferr, float * berr, float * work, 
                       int * iwork, int * info);

extern "C" void PREFIX spotrf_(Epetra_fcd, int* n, float* a, int* lda, int* info); 
extern "C" void PREFIX spotrs_(Epetra_fcd, int * n, int * nrhs, float* a,
                       int* lda, float*x , int* ldx, int* info);
extern "C" void PREFIX spotri_(Epetra_fcd, int* n, float* a, int* lda, int* info); 
extern "C" void PREFIX spocon_(Epetra_fcd, int* n, float* a, int* lda, 
                       float * anorm, float * rcond, float * work,
                       int * iwork, int* info); 
extern "C" void PREFIX sposv_(Epetra_fcd, int * n, int * nrhs, float* a,
                       int* lda, float*x , int* ldx, int* info);
extern "C" void PREFIX spoequ_(int* n, float* a, int* lda, float * s, float * scond,
                       float * amax, int* info); 

extern "C" void PREFIX sporfs_(Epetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, 
                       float * b, int * ldb, float * x, int * ldx, 
                       float * ferr, float * berr, float * work, 
                       int * iwork, int * info);

extern "C" void PREFIX sposvx_(Epetra_fcd, Epetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, Epetra_fcd, 
                       float * s, float * b, int * ldb, float * x, int * ldx, 
                       float * rcond, float * ferr, float * berr, float * work, 
                       int * iwork, int * info);

// Double precision LAPACK eigen solvers
extern "C" void PREFIX dgels_(Epetra_fcd ch, int*, int*, int*,
                       double*, int*, double*, int*, double*, int*, int*);

extern "C" void PREFIX dgeev_(Epetra_fcd, Epetra_fcd, int*, double*, int*,
                      double*, double*, double*, int*, double*, int*,
                      double*, int*, int*);
           
extern "C" void PREFIX  dgehrd_(int * n, int * ilo, int * ihi, double * A,
                        int * lda, double * tau, double * work, int * lwork,
                        int * info);

extern "C" void PREFIX  dhseqr_(Epetra_fcd job, Epetra_fcd, int * n, int * ilo, int * ihi,
                        double * h, int * ldh, double * wr, double * wi, double * z,
                        int * ldz, double * work, int * lwork, int * info);

extern "C" void PREFIX  dorghr_(int * n, int * ilo, int * ihi, double * a, int * lda, double * tau,
                        double * work, int * lwork, int * info);
                        
extern "C" void PREFIX  dormhr_(Epetra_fcd, Epetra_fcd, int * m, int * n, int * ilo,
                        int * ihi, double * a, int * lda, double * tau, double * c,
                        int * ldc, double * work, int * lwork, int * info);

extern "C" void PREFIX  dtrevc_(Epetra_fcd, Epetra_fcd, int * select, int * n, double * t,
                        int * ldt, double *vl, int * ldvl, double * vr, int * ldvr,
                        int * mm, int * m, double * work, int * info); 

extern "C" void PREFIX  dtrexc_(Epetra_fcd, int * n, double * t, int * ldt, double * q,
                        int * ldq, int * ifst, int * ilst, double * work, int * info);

extern "C" double PREFIX  dlamch_(Epetra_fcd);


// Single precision LAPACK eigen solvers
extern "C" void PREFIX sgels_(Epetra_fcd, int*, int*, int*,
                       float*, int*, float*, int*, float*, int*, int*);

extern "C" void PREFIX sgeev_(Epetra_fcd, Epetra_fcd, int*, float*, int*,
                      float*, float*, float*, int*, float*, int*,
                      float*, int*, int*);

extern "C" void PREFIX  sgehrd_(int * n, int * ilo, int * ihi, float * A,
                        int * lda, float * tau, float * work, int * lwork,
                        int * info);

extern "C" void PREFIX  shseqr_(Epetra_fcd job, Epetra_fcd, int * n, int * ilo, int * ihi,
                        float * h, int * ldh, float * wr, float * wi, float * z,
                        int * ldz, float * work, int * lwork, int * info);

extern "C" void PREFIX  sorghr_(int * n, int * ilo, int * ihi, float * a, int * lda, float * tau,
                        float * work, int * lwork, int * info);
                        
extern "C" void PREFIX  sormhr_(Epetra_fcd, Epetra_fcd, int * m, int * n, int * ilo,
                        int * ihi, float * a, int * lda, float * tau, float * c,
                        int * ldc, float * work, int * lwork, int * info);

extern "C" void PREFIX  strevc_(Epetra_fcd, Epetra_fcd, int * select, int * n, float * t,
                        int * ldt, float *vl, int * ldvl, float * vr, int * ldvr,
                        int * mm, int * m, float * work, int * info); 

extern "C" void PREFIX  strexc_(Epetra_fcd, int * n, float * t, int * ldt, float * q,
                        int * ldq, int * ifst, int * ilst, float * work, int * info);

extern "C" float PREFIX  slamch_(Epetra_fcd);
