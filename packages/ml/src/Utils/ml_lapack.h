/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_BLAS_LAPACK_WRAPPERS_H
#define ML_BLAS_LAPACK_WRAPPERS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/* Wrappers for blas & lapack routines.   This file is the result of merging
 * Epetra_BLAS_wrappers.h & Epetra_LAPACK_wrappers.h and a few cosmetic
 * changes.  */

#include "ml_common.h"

#ifdef f2c_i2
/* for -i2 */
typedef short ftnlen;
#else
typedef long ftnlen;
#endif

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_BLAS_MANGLE is defined undefine it because we want to redefine */

#if defined(CRAY_T3X)

#include "fortran.h"
#define ml_fcd fcd
#define PREFIX

#define DGBSV_F77   F77_BLAS_MANGLE(sgbsv,SGBSV)
#define DGTTRF_F77  F77_BLAS_MANGLE(sgttrf,SGTTRF)
#define DGTTS2_F77  F77_BLAS_MANGLE(sgtts2,SGTTS2)
#define DGETRF_F77  F77_BLAS_MANGLE(sgetrf,SGETRF)
#define DGETRS_F77  F77_BLAS_MANGLE(sgetrs,SGETRS)
#define DGESVD_F77  F77_BLAS_MANGLE(sgesvd,SGESVD)
#define DPOTRF_F77  F77_BLAS_MANGLE(spotrf,SPOTRF)
#define DPOTRS_F77  F77_BLAS_MANGLE(spotrs,SPOTRS)
#define DPOTRI_F77  F77_BLAS_MANGLE(spotri,SPOTRI)
#define DPOCON_F77  F77_BLAS_MANGLE(spocon,SPOCON)
#define DPOSV_F77   F77_BLAS_MANGLE(sposv,SPOSV)
#define DPOEQU_F77  F77_BLAS_MANGLE(spoequ,SPOEQU)
#define DPORFS_F77  F77_BLAS_MANGLE(sporfs,SPORFS)
#define DPOSVX_F77  F77_BLAS_MANGLE(sposvx,SPOSVX)
#define DGELS_F77   F77_BLAS_MANGLE(sgels,SGELS)
#define DGEEV_F77   F77_BLAS_MANGLE(sgeev,SGEEV)
#define DGEHRD_F77  F77_BLAS_MANGLE(sgehrs,SGEHRS)
#define DHSEQR_F77  F77_BLAS_MANGLE(shseqr,SHSEQR)
#define DORGHR_F77  F77_BLAS_MANGLE(sorghr,SORGHR)
#define DORMHR_F77  F77_BLAS_MANGLE(sormhr,SORMHR)
#define DTREVC_F77  F77_BLAS_MANGLE(strevc,STREVC)
#define DTREXC_F77  F77_BLAS_MANGLE(strexc,STREXC)
#define DGELSS_F77  F77_BLAS_MANGLE(sgelss,SGELSS)
#define DSTEV_F77   F77_BLAS_MANGLE(sstev,SSTEV)
#define DGEQPF_F77  F77_BLAS_MANGLE(sgeqpf,SGEQPF)
#define DGEQRF_F77  F77_BLAS_MANGLE(sgeqrf,SGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(sorgqr,SORGQR)
#define DSTEQR_F77  F77_BLAS_MANGLE(ssteqr,SSTEQR)

#endif
#if defined(INTEL_CXML)

#define ml_fcd char *, unsigned int
#define PREFIX __stdcall

#define DGEEV_F77   F77_BLAS_MANGLE(dgeev,DGEEV)
#define DGEHRD_F77  F77_BLAS_MANGLE(dgehrd,DGEHRD)
#define DGELS_F77   F77_BLAS_MANGLE(dgels,DGELS)
#define DGBSV_F77   F77_BLAS_MANGLE(dgbsv,DGBSV)
#define DGTTRF_F77  F77_BLAS_MANGLE(dgttrf,DGTTRF)
#define DGTTS2_F77  F77_BLAS_MANGLE(dgtts2,DGTTS2)
#define DGETRF_F77  F77_BLAS_MANGLE(dgetrf,DGETRF)
#define DGETRS_F77  F77_BLAS_MANGLE(dgetrs,DGETRS)
#define DGESVD_F77  F77_BLAS_MANGLE(dgesvd,DGESVD)
#define DHSEQR_F77  F77_BLAS_MANGLE(dhseqr,DHSEQR)
#define DORGHR_F77  F77_BLAS_MANGLE(dorghr,DORGHR)
#define DORMHR_F77  F77_BLAS_MANGLE(dormhr,DORMHR)
#define DPOTRF_F77  F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define DPOTRS_F77  F77_BLAS_MANGLE(dpotrs,DPOTRS)
#define DPOTRI_F77  F77_BLAS_MANGLE(dpotri,DPOTRI)
#define DPOCON_F77  F77_BLAS_MANGLE(dpocon,DPOCON)
#define DPOSV_F77   F77_BLAS_MANGLE(dposv,DPOSV)
#define DPOEQU_F77  F77_BLAS_MANGLE(dpoequ,DPOEQU)
#define DPORFS_F77  F77_BLAS_MANGLE(dporfs,DPORFS)
#define DPOSVX_F77  F77_BLAS_MANGLE(dposvx,DPOSVX)
#define DTREVC_F77  F77_BLAS_MANGLE(dtrevc,DTREVC)
#define DTREXC_F77  F77_BLAS_MANGLE(dtrexc,DTREXC)
#define DGELSS_F77  F77_BLAS_MANGLE(dgelss,DGELSS)
#define DSTEV_F77   F77_BLAS_MANGLE(dstev,DSTEV)
#define DGEQPF_F77  F77_BLAS_MANGLE(dgeqpf,DGEQPF)
#define DGEQRF_F77  F77_BLAS_MANGLE(dgeqrf,DGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(dorgqr,DORGQR)
#define DSTEQR_F77  F77_BLAS_MANGLE(dsteqr,DSTEQR)

#endif
#if defined(INTEL_MKL)

#define ml_fcd char *
#define PREFIX

#define DGBSV_F77   F77_BLAS_MANGLE(dgbsv,DGBSV)
#define DGTTRF_F77  F77_BLAS_MANGLE(dgttrf,DGTTRF)
#define DGTTS2_F77  F77_BLAS_MANGLE(dgtts2,DGTTS2)
#define DGETRF_F77  F77_BLAS_MANGLE(dgetrf,DGETRF)
#define DGETRS_F77  F77_BLAS_MANGLE(dgetrs,DGETRS)
#define DGESVD_F77  F77_BLAS_MANGLE(dgesvd,DGESVD)
#define DPOTRF_F77  F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define DPOTRS_F77  F77_BLAS_MANGLE(dpotrs,DPOTRS)
#define DPOTRI_F77  F77_BLAS_MANGLE(dpotri,DPOTRI)
#define DPOCON_F77  F77_BLAS_MANGLE(dpocon,DPOCON)
#define DPOSV_F77   F77_BLAS_MANGLE(dposv,DPOSV)
#define DPOEQU_F77  F77_BLAS_MANGLE(dpoequ,DPOEQU)
#define DPORFS_F77  F77_BLAS_MANGLE(dporfs,DPORFS)
#define DPOSVX_F77  F77_BLAS_MANGLE(dposvx,DPOSVX)
#define DGELS_F77   F77_BLAS_MANGLE(dgels,DGELS)
#define DGEEV_F77   F77_BLAS_MANGLE(dgeev,DGEEV)
#define DGEHRD_F77  F77_BLAS_MANGLE(dgehrd,DGEHRD)
#define DHSEQR_F77  F77_BLAS_MANGLE(dhseqr,DHSEQR)
#define DORGHR_F77  F77_BLAS_MANGLE(dorghr,DORGHR)
#define DORMHR_F77  F77_BLAS_MANGLE(dormhr,DORMHR)
#define DTREVC_F77  F77_BLAS_MANGLE(dtrevc,DTREVC)
#define DTREXC_F77  F77_BLAS_MANGLE(dtrexc,DTREXC)
#define DGELSS_F77  F77_BLAS_MANGLE(dgelss,DGELSS)
#define DSTEV_F77   F77_BLAS_MANGLE(dstev,DSTEV)
#define DGEQPF_F77  F77_BLAS_MANGLE(dgeqpf,DGEQPF)
#define DGEQRF_F77  F77_BLAS_MANGLE(dgeqrf,DGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(dorgqr,DORGQR)
#define DSTEQR_F77  F77_BLAS_MANGLE(dsteqr,DSTEQR)

#endif

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif

#define F77_BLAS_MANGLE(lcase,UCASE) UCASE

#else

#define ml_fcd char *
#define PREFIX

/* Use autoconf's definition of F77_BLAS_MANGLE
   unless using old make system */

#ifdef TRILINOS_NO_CONFIG_H

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_BLAS_MANGLE(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_BLAS_MANGLE(lcase,UCASE) lcase ## _
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */
#endif /* TRILINOS_NO_CONFIG_H */

#define DGBSV_F77   F77_BLAS_MANGLE(dgbsv,DGBSV)
#define DGTTRF_F77  F77_BLAS_MANGLE(dgttrf,DGTTRF)
#define DGTTS2_F77  F77_BLAS_MANGLE(dgtts2,DGTTS2)
#define DGETRF_F77  F77_BLAS_MANGLE(dgetrf,DGETRF)
#define DGETRS_F77  F77_BLAS_MANGLE(dgetrs,DGETRS)
#define DGESVD_F77  F77_BLAS_MANGLE(dgesvd,DGESVD)
#define DPOTRF_F77  F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define DPOTRS_F77  F77_BLAS_MANGLE(dpotrs,DPOTRS)
#define DPOTRI_F77  F77_BLAS_MANGLE(dpotri,DPOTRI)
#define DPOCON_F77  F77_BLAS_MANGLE(dpocon,DPOCON)
#define DPOSV_F77   F77_BLAS_MANGLE(dposv,DPOSV)
#define DPOEQU_F77  F77_BLAS_MANGLE(dpoequ,DPOEQU)
#define DPORFS_F77  F77_BLAS_MANGLE(dporfs,DPORFS)
#define DPOSVX_F77  F77_BLAS_MANGLE(dposvx,DPOSVX)
#define DGELS_F77   F77_BLAS_MANGLE(dgels,DGELS)
#define DGEEV_F77   F77_BLAS_MANGLE(dgeev,DGEEV)
#define DGEHRD_F77  F77_BLAS_MANGLE(dgehrd,DGEHRD)
#define DHSEQR_F77  F77_BLAS_MANGLE(dhseqr,DHSEQR)
#define DORGHR_F77  F77_BLAS_MANGLE(dorghr,DORGHR)
#define DORMHR_F77  F77_BLAS_MANGLE(dormhr,DORMHR)
#define DTREVC_F77  F77_BLAS_MANGLE(dtrevc,DTREVC)
#define DTREXC_F77  F77_BLAS_MANGLE(dtrexc,DTREXC)
#define DGELSS_F77  F77_BLAS_MANGLE(dgelss,DGELSS)
#define DSTEV_F77   F77_BLAS_MANGLE(dstev,DSTEV)
#define DGEQPF_F77  F77_BLAS_MANGLE(dgeqpf,DGEQPF)
#define DGEQRF_F77  F77_BLAS_MANGLE(dgeqrf,DGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(dorgqr,DORGQR)
#define DSTEQR_F77  F77_BLAS_MANGLE(dsteqr,DSTEQR)

#endif

#define DGETRI_F77  F77_BLAS_MANGLE(dgetri,DGETRI)
#define DGERFS_F77  F77_BLAS_MANGLE(dgerfs,DGERFS)
#define DGECON_F77  F77_BLAS_MANGLE(dgecon,DGECON)
#define DGESVX_F77  F77_BLAS_MANGLE(dgesvx,DGESVX)
#define DGESV_F77   F77_BLAS_MANGLE(dgesv,DGESV)
#define DGEEQU_F77  F77_BLAS_MANGLE(dgeequ,DGEEQU)
#define DLAMCH_F77  F77_BLAS_MANGLE(dlamch,DLAMCH)
#define DGELSS_F77  F77_BLAS_MANGLE(dgelss,DGELSS)
#define DSTEV_F77   F77_BLAS_MANGLE(dstev,DSTEV)
#define DGEQPF_F77  F77_BLAS_MANGLE(dgeqpf,DGEQPF)
#define DGEQRF_F77  F77_BLAS_MANGLE(dgeqrf,DGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(dorgqr,DORGQR)
#define DSTEQR_F77  F77_BLAS_MANGLE(dsteqr,DSTEQR)

#define SGBSV_F77   F77_BLAS_MANGLE(sgbsv,SGBSV)
#define SGTTRF_F77  F77_BLAS_MANGLE(sgttrf,SGTTRF)
#define SGTTS2_F77  F77_BLAS_MANGLE(sgtts2,SGTTS2)
#define SGETRF_F77  F77_BLAS_MANGLE(sgetrf,SGETRF)
#define SGETRS_F77  F77_BLAS_MANGLE(sgetrs,SGETRS)
#define SGESVD_F77  F77_BLAS_MANGLE(sgesvd,SGESVD)
#define SGETRI_F77  F77_BLAS_MANGLE(sgetri,SGETRI)
#define SGERFS_F77  F77_BLAS_MANGLE(sgerfs,SGERFS)
#define SGECON_F77  F77_BLAS_MANGLE(sgecon,SGECON)
#define SGESVX_F77  F77_BLAS_MANGLE(sgesvx,SGESVX)
#define SGESV_F77   F77_BLAS_MANGLE(sgesv,SGESV)
#define SGEEQU_F77  F77_BLAS_MANGLE(sgeequ,SGEEQU)
#define SPOTRF_F77  F77_BLAS_MANGLE(spotrf,SPOTRF)
#define SPOTRS_F77  F77_BLAS_MANGLE(spotrs,SPOTRS)
#define SPOTRI_F77  F77_BLAS_MANGLE(spotri,SPOTRI)
#define SPOCON_F77  F77_BLAS_MANGLE(spocon,SPOCON)
#define SPOSV_F77   F77_BLAS_MANGLE(sposv,SPOSV)
#define SPOEQU_F77  F77_BLAS_MANGLE(spoequ,SPOEQU)
#define SPORFS_F77  F77_BLAS_MANGLE(sporfs,SPORFS)
#define SPOSVX_F77  F77_BLAS_MANGLE(sposvx,SPOSVX)
#define SGELS_F77   F77_BLAS_MANGLE(sgels,SGELS)
#define SGEEV_F77   F77_BLAS_MANGLE(sgeev,SGEEV)
#define SGEHRD_F77  F77_BLAS_MANGLE(sgehrd,SGEHRD)
#define SHSEQR_F77  F77_BLAS_MANGLE(shseqr,SHSEQR)
#define SORGHR_F77  F77_BLAS_MANGLE(sorghr,SORGHR)
#define SORMHR_F77  F77_BLAS_MANGLE(sormhr,SORMHR)
#define STREVC_F77  F77_BLAS_MANGLE(strevc,STREVC)
#define STREXC_F77  F77_BLAS_MANGLE(strexc,STREXC)
#define SLAMCH_F77  F77_BLAS_MANGLE(slamch,SLAMCH)
#define SGELSS_F77  F77_BLAS_MANGLE(sgelss,SGELSS)
#define SSTEV_F77   F77_BLAS_MANGLE(sstev,SSTEV)
#define SGEQPF_F77  F77_BLAS_MANGLE(sgeqpf,SGEQPF)

#define DLASWP_F77  F77_BLAS_MANGLE(dlaswp,DLASWP)
#define DLAIC1_F77  F77_BLAS_MANGLE(dlaic1,DLAIC1)


#ifdef __cplusplus
extern "C" {
#endif

/* Aztec-2.1 already defines these variables in header files */

#ifndef HAVE_ML_AZTEC2_1

  /* Double precision LAPACK linear solvers */
void PREFIX DGBSV_F77(int *n, int *m, int *i, int *k, double *a, int *j, int *ipiv, double * b, int *nn, int *info);
void PREFIX DGETRF_F77(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
void PREFIX DGTTRF_F77(int* leng, double* a, double *b, double *c, double *d, int *e, int *info);
void PREFIX DGTTS2_F77(ml_fcd, int* n, int* nrhs, double* a, double* b,
                       double* c,double* d,int* e,double* x, int* ldx);
void PREFIX DGETRS_F77(ml_fcd, int* n, int* nrhs, double* a,
                       int* lda, int* ipiv, double* x , int* ldx, int* info);
void PREFIX DGESVD_F77(ml_fcd, ml_fcd, int* m, int* n, double* a,
                       int* lda, double* s, double* u, int* ldu, double* vt, int* ldut,
		       double* work, int* lwork, int* info);
void PREFIX DGETRI_F77(int* n, double* a, int* lda, int*ipiv, double * work , int* lwork, int* info);
void PREFIX DGECON_F77(ml_fcd norm, int* n, double* a, int* lda,
                       double *anorm, double * rcond, double * work,
                       int * iwork, int* info);
void PREFIX DGESV_F77(int * n, int * nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void PREFIX DGEEQU_F77(int* m, int* n, double* a, int* lda, double * r, double * c,
			double * rowcnd, double * colcnd,
                       double * amax, int* info);
void PREFIX DGERFS_F77(ml_fcd, int * n, int * nrhs, double * a,
                       int * lda, double * af, int * ldaf, int*ipiv,
                       double * b, int * ldb, double * x, int * ldx,
                       double * ferr, double * berr, double * work,
                       int * iwork, int * info);
void PREFIX DGESVX_F77(ml_fcd, ml_fcd, int * n, int * nrhs, double * a,
                       int * lda, double * af, int * ldaf, int*ipiv, ml_fcd,
                       double * r, double *c, double * b, int * ldb, double * x, int * ldx,
                       double * rcond, double * ferr, double * berr, double * work,
                       int * iwork, int * info);

void PREFIX DPOTRF_F77(ml_fcd, int* n, double* a, int* lda, int* info);
void PREFIX DPOTRS_F77(ml_fcd, int * n, int * nrhs, double* a,
                       int* lda, double*x , int* ldx, int* info);
void PREFIX DPOTRI_F77(ml_fcd, int* n, double* a, int* lda, int* info);
void PREFIX DPOCON_F77(ml_fcd, int* n, double* a, int* lda,
                       double * anorm, double * rcond, double * work,
                       int * iwork, int* info);
void PREFIX DPOSV_F77(ml_fcd, int * n, int * nrhs, double* a,
                       int* lda, double*x , int* ldx, int* info);
void PREFIX DPOEQU_F77(int* n, double* a, int* lda, double * s, double * scond,
                       double * amax, int* info);

void PREFIX DPORFS_F77(ml_fcd, int * n, int * nrhs, double * a,
                       int * lda, double * af, int * ldaf,
                       double * b, int * ldb, double * x, int * ldx,
                       double * ferr, double * berr, double * work,
                       int * iwork, int * info);

void PREFIX DPOSVX_F77(ml_fcd, ml_fcd, int * n, int * nrhs, double * a,
                       int * lda, double * af, int * ldaf, ml_fcd,
                       double * s, double * b, int * ldb, double * x, int * ldx,
                       double * rcond, double * ferr, double * berr, double * work,
                       int * iwork, int * info);

void PREFIX DGELSS_F77(int * m, int * n, int * nrhs, double * a,
				   int * lda, double * b, int * ldb, double * s,
				   double * rcond, int * rank, double * work,
				   int * lwork, int * info);

void PREFIX DGEQPF_F77(int * m, int * n, double * a,
				   int * lda, int * jpvt, double * tau,
				   double * work, int * info);
/*
#define DGEQRF_F77  F77_BLAS_MANGLE(sgeqrf,SGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(sorgqr,SORGQR)
*/
void PREFIX DGEQRF_F77(int *, int *, double *, int *,
                 double *, double *, int *, int *);

void PREFIX DORGQR_F77(int *m, int *n, int *k, double * a,
                 int *lda, double *tau, double *work, int *lwork,
                 int *info);
void PREFIX DSTEQR_F77(ml_fcd, int *, double *, double *, double *, int *,
                 double *, int *);

  /* Single precision LAPACK linear solvers*/
void PREFIX SGETRF_F77(int* m, int* n, float* a, int* lda, int* ipiv, int* info);
void PREFIX SGETRS_F77(ml_fcd, int* n, int* nrhs, float* a,
                       int* lda, int* ipiv, float* x , int* ldx, int* info);
void PREFIX SGESVD_F77(ml_fcd, ml_fcd, int* m, int* n, float* a,
		       int* lda, float* s, float* u, int* ldu, float* vt, int* ldut,
		       float* work, int* lwork, int* info);
void PREFIX SGETRI_F77(int* n, float* a, int* lda, int*ipiv, float * work , int* lwork, int* info);
void PREFIX SGECON_F77(ml_fcd norm, int* n, float* a, int* lda,
                       float * anorm, float * rcond, float * work,
                       int * iwork, int* info);
void PREFIX SGESV_F77(int * n, int * nrhs, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
void PREFIX SGEEQU_F77(int* m, int* n, float* a, int* lda, float * r, float * c,
			float * rowcnd, float * colcnd,
			float * amax, int* info);
void PREFIX SGERFS_F77(ml_fcd, int * n, int * nrhs, float * a,
                       int * lda, float * af, int * ldaf, int*ipiv,
                       float * b, int * ldb, float * x, int * ldx,
                       float * ferr, float * berr, float * work,
                       int * iwork, int * info);
void PREFIX SGESVX_F77(ml_fcd, ml_fcd, int * n, int * nrhs, float * a,
                       int * lda, float * af, int * ldaf, int*ipiv, ml_fcd,
                       float * r, float *c, float * b, int * ldb, float * x, int * ldx,
                       float * rcond, float * ferr, float * berr, float * work,
                       int * iwork, int * info);

void PREFIX SPOTRF_F77(ml_fcd, int* n, float* a, int* lda, int* info);
void PREFIX SPOTRS_F77(ml_fcd, int * n, int * nrhs, float* a,
                       int* lda, float*x , int* ldx, int* info);
void PREFIX SPOTRI_F77(ml_fcd, int* n, float* a, int* lda, int* info);
void PREFIX SPOCON_F77(ml_fcd, int* n, float* a, int* lda,
                       float * anorm, float * rcond, float * work,
                       int * iwork, int* info);
void PREFIX SPOSV_F77(ml_fcd, int * n, int * nrhs, float* a,
                       int* lda, float*x , int* ldx, int* info);
void PREFIX SPOEQU_F77(int* n, float* a, int* lda, float * s, float * scond,
                       float * amax, int* info);

void PREFIX SPORFS_F77(ml_fcd, int * n, int * nrhs, float * a,
                       int * lda, float * af, int * ldaf,
                       float * b, int * ldb, float * x, int * ldx,
                       float * ferr, float * berr, float * work,
                       int * iwork, int * info);

void PREFIX SPOSVX_F77(ml_fcd, ml_fcd, int * n, int * nrhs, float * a,
                       int * lda, float * af, int * ldaf, ml_fcd,
                       float * s, float * b, int * ldb, float * x, int * ldx,
                       float * rcond, float * ferr, float * berr, float * work,
                       int * iwork, int * info);

void PREFIX SGELSS_F77(int * m, int * n, int * nrhs, float * a,
                       int * lda, float * b, int * ldb, float * s,
                       float * rcond, int * rank, float * work,
                       int * lwork, int * info);

void PREFIX SGEQPF_F77(int * m, int * n, float * a,
                       int * lda, int * jpvt, float * tau,
                       float * work, int * info);

  /* Double precision LAPACK eigen solvers*/
void PREFIX DGELS_F77(ml_fcd ch, int*, int*, int*,
                       double*, int*, double*, int*, double*, int*, int*);

void PREFIX DGEEV_F77(ml_fcd, ml_fcd, int*, double*, int*,
                      double*, double*, double*, int*, double*, int*,
                      double*, int*, int*);

void PREFIX DGEHRD_F77(int * n, int * ilo, int * ihi, double * A,
                        int * lda, double * tau, double * work, int * lwork,
                        int * info);

void PREFIX DHSEQR_F77(ml_fcd job, ml_fcd, int * n, int * ilo, int * ihi,
                        double * h, int * ldh, double * wr, double * wi, double * z,
                        int * ldz, double * work, int * lwork, int * info);

void PREFIX DORGHR_F77(int * n, int * ilo, int * ihi, double * a, int * lda, double * tau,
                        double * work, int * lwork, int * info);

void PREFIX DORMHR_F77(ml_fcd, ml_fcd, int * m, int * n, int * ilo,
                        int * ihi, double * a, int * lda, double * tau, double * c,
                        int * ldc, double * work, int * lwork, int * info);

void PREFIX DTREVC_F77(ml_fcd, ml_fcd, int * select, int * n, double * t,
                        int * ldt, double *vl, int * ldvl, double * vr, int * ldvr,
                        int * mm, int * m, double * work, int * info);

void PREFIX DTREXC_F77(ml_fcd, int * n, double * t, int * ldt, double * q,
                        int * ldq, int * ifst, int * ilst, double * work, int * info);

void PREFIX DSTEV_F77(ml_fcd jobz, int * n, double * d,
                      double * e, double * z, int * ldz,
                      double * work, int * info);

double PREFIX DLAMCH_F77(ml_fcd);


  /* Single precision LAPACK eigen solvers*/
void PREFIX SGELS_F77(ml_fcd, int*, int*, int*,
                       float*, int*, float*, int*, float*, int*, int*);

void PREFIX SGEEV_F77(ml_fcd, ml_fcd, int*, float*, int*,
                      float*, float*, float*, int*, float*, int*,
                      float*, int*, int*);

void PREFIX SGEHRD_F77(int * n, int * ilo, int * ihi, float * A,
                        int * lda, float * tau, float * work, int * lwork,
                        int * info);

void PREFIX SHSEQR_F77(ml_fcd job, ml_fcd, int * n, int * ilo, int * ihi,
                        float * h, int * ldh, float * wr, float * wi, float * z,
                        int * ldz, float * work, int * lwork, int * info);

void PREFIX SORGHR_F77(int * n, int * ilo, int * ihi, float * a, int * lda, float * tau,
                        float * work, int * lwork, int * info);

void PREFIX SORMHR_F77(ml_fcd, ml_fcd, int * m, int * n, int * ilo,
                        int * ihi, float * a, int * lda, float * tau, float * c,
                        int * ldc, float * work, int * lwork, int * info);

void PREFIX STREVC_F77(ml_fcd, ml_fcd, int * select, int * n, float * t,
                        int * ldt, float *vl, int * ldvl, float * vr, int * ldvr,
                        int * mm, int * m, float * work, int * info);

void  PREFIX STREXC_F77(ml_fcd, int * n, float * t, int * ldt, float * q,
                        int * ldq, int * ifst, int * ilst, float * work, int * info);

void PREFIX SSTEV_F77(ml_fcd jobz, int * n, float * d,
                      float * e, float * z, int * ldz,
                      float * work, int * info);

float  PREFIX SLAMCH_F77(ml_fcd);

#endif

#ifdef __cplusplus
}
#endif


/* Define fcd (Fortran ml_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)


#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define ml_fcd fcd

#define DASUM_F77   F77_BLAS_MANGLE(sasum,SASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(saxpy,SAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(scopy,SCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(sdot,SDOT)
#define DNRM2_F77   F77_BLAS_MANGLE(snrm2,SNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(sscal,SSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(isamax,ISAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(sgemv,SGEMV)
#define DGER_F77    F77_BLAS_MANGLE(sger,SGER)
#define DTRMV_F77   F77_BLAS_MANGLE(strmv,STRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(sgemm,SGEMM)
#define DSYMM_F77   F77_BLAS_MANGLE(ssymm,SSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(strmm,STRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(strsm,STRSM)

#elif defined(INTEL_CXML)

#define PREFIX __stdcall
#define ml_fcd char *, unsigned int

#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)
#define DNRM2_F77   F77_BLAS_MANGLE(dnrm2,DNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGER_F77    F77_BLAS_MANGLE(dger,DGER)
#define DTRMV_F77   F77_BLAS_MANGLE(dtrmv,DTRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM)
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)


#elif defined(INTEL_MKL)

#define PREFIX
#define ml_fcd char *

#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)
#define DNRM2_F77   F77_BLAS_MANGLE(dnrm2,DNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGER_F77    F77_BLAS_MANGLE(dger,DGER)
#define DTRMV_F77   F77_BLAS_MANGLE(dtrmv,DTRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM)
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)


#endif

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_BLAS_MANGLE is defined undefine it because we want to redefine */

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif


#define F77_BLAS_MANGLE(lcase,UCASE) UCASE

#else /* Define ml_fcd for all other machines */

#define PREFIX
#define ml_fcd char *

/* Use autoconf's definition of F77_BLAS_MANGLE
   unless using old make system */

#ifdef TRILINOS_NO_CONFIG_H

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_BLAS_MANGLE(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_BLAS_MANGLE(lcase,UCASE) lcase ## _
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */
#endif /* TRILINOS_NO_CONFIG_H */

#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)
#define DNRM2_F77   F77_BLAS_MANGLE(dnrm2,DNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGER_F77    F77_BLAS_MANGLE(dger,DGER)
#define DTRMV_F77   F77_BLAS_MANGLE(dtrmv,DTRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM)
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)


#endif


#define SSCAL_F77   F77_BLAS_MANGLE(sscal,SSCAL)
#define SCOPY_F77   F77_BLAS_MANGLE(scopy,SCOPY)
#define SAXPY_F77   F77_BLAS_MANGLE(saxpy,SAXPY)
#define SDOT_F77    F77_BLAS_MANGLE(sdot,SDOT)
#define SNRM2_F77   F77_BLAS_MANGLE(snrm2,SNRM2)
#define SASUM_F77   F77_BLAS_MANGLE(sasum,SASUM)
#define ISAMAX_F77  F77_BLAS_MANGLE(isamax,ISAMAX)

#define SGEMV_F77   F77_BLAS_MANGLE(sgemv,SGEMV)
#define SGER_F77    F77_BLAS_MANGLE(sger,SGER)
#define STRMV_F77   F77_BLAS_MANGLE(strmv,STRMV)
#define SGEMM_F77   F77_BLAS_MANGLE(sgemm,SGEMM)
#define SSYMM_F77   F77_BLAS_MANGLE(ssymm,SSYMM)
#define STRMM_F77   F77_BLAS_MANGLE(strmm,STRMM)
#define STRSM_F77   F77_BLAS_MANGLE(strsm,STRSM)

/* Explicitly define each F77 name for all BLAS kernels */

#ifdef __cplusplus
extern "C" {
#endif

/* Aztec-2.1 already contains these definitions in the header files */
#ifndef HAVE_ML_AZTEC2_1
/* Double precision BLAS 1 */
double PREFIX DASUM_F77(int* n, double x[], int* incx);
void PREFIX DAXPY_F77(int* n, double* alpha, double x[], int* incx, double y[], int* incy);
void PREFIX DCOPY_F77(int* n, double *x, int* incx, double *y, int* incy);
double PREFIX DDOT_F77(int* n, double x[], int* incx, double y[], int* incy);
double PREFIX DNRM2_F77(int* n, double x[], int* incx);
void PREFIX DSCAL_F77(int* n, double* alpha, double *x, int* incx);
int PREFIX IDAMAX_F77(int* n, double *x, int* incx);

/* Single precision BLAS 1 */
float PREFIX SASUM_F77(int* n, float x[], int* incx);
void PREFIX SAXPY_F77(int* n, float* alpha, float x[], int* incx, float y[], int* incy);
void PREFIX SCOPY_F77(int* n, float *x, int* incx, float *y, int* incy);
float PREFIX SDOT_F77(int* n, float x[], int* incx, float y[], int* incy);
float PREFIX SNRM2_F77(int* n, float x[], int* incx);
void PREFIX SSCAL_F77(int* n, float* alpha, float *x, int* incx);
int PREFIX ISAMAX_F77(int* n, float *x, int* incx);

/* Double precision BLAS 2 */
void PREFIX DGEMV_F77(ml_fcd, int* m, int* n, double* alpha, double A[], int* lda,
		       double x[], int* incx, double* beta, double y[], int* incy);
void PREFIX DTRMV_F77(ml_fcd, ml_fcd, ml_fcd, int *n,
		      double *a, int *lda, double *x, int *incx);
void PREFIX DGER_F77(int *m, int *n, double *alpha, double *x, int *incx, double *y,
		     int *incy, double *a, int *lda);


/* Single precision BLAS 2 */
void PREFIX SGEMV_F77(ml_fcd, int* m, int* n, float* alpha, float A[], int* lda,
		       float x[], int* incx, float* beta, float y[], int* incy);
void PREFIX STRMV_F77(ml_fcd, ml_fcd, ml_fcd, int *n,
		      float *a, int *lda, float *x, int *incx);
void PREFIX SGER_F77(int *m, int *n, float *alpha, float *x, int *incx, float *y,
		     int *incy, float *a, int *lda);

/* Double precision BLAS 3 */
void PREFIX DGEMM_F77(ml_fcd, ml_fcd, int *m, int *
		      n, int *k, double *alpha, double *a, int *lda,
		      double *b, int *ldb, double *beta, double *c, int *ldc);
void PREFIX DSYMM_F77(ml_fcd, ml_fcd, int *m, int * n,
		      double *alpha, double *a, int *lda,
		      double *b, int *ldb, double *beta, double *c, int *ldc);
void PREFIX DTRMM_F77(ml_fcd, ml_fcd, ml_fcd, ml_fcd,
		      int *m, int *n, double *alpha, double *a, int * lda, double *b, int *ldb);
void PREFIX DTRSM_F77(ml_fcd, ml_fcd, ml_fcd, ml_fcd,
		      int *m, int *n, double *alpha, double *a, int *
		      lda, double *b, int *ldb);

/* Single precision BLAS 3 */
void PREFIX SGEMM_F77(ml_fcd, ml_fcd, int *m, int *
		      n, int *k, float *alpha, float *a, int *lda,
		      float *b, int *ldb, float *beta, float *c, int *ldc);
void PREFIX SSYMM_F77(ml_fcd, ml_fcd, int *m, int * n,
		      float *alpha, float *a, int *lda,
		      float *b, int *ldb, float *beta, float *c, int *ldc);
void PREFIX STRMM_F77(ml_fcd, ml_fcd, ml_fcd, ml_fcd,
		      int *m, int *n, float *alpha, float *a, int * lda, float *b, int *ldb);
void PREFIX STRSM_F77(ml_fcd, ml_fcd, ml_fcd, ml_fcd,
		      int *m, int *n, float *alpha, float *a, int *
		      lda, float *b, int *ldb);

void PREFIX XERBLA_F77(ml_fcd, int *info);

#endif

#ifdef __cplusplus
}
#endif

#endif /* ML_BLAS_LAPACK_WRAPPERS_H */
