// 27-May-2002 General cleanup. Changed method names to fit namingConvention (nothing changed).

#ifndef _TPETRA_LAPACK_WRAPPERS_H_
#define _TPETRA_LAPACK_WRAPPERS_H_

/*! \file Tpetra_LAPACK_wrappers.h
    \brief It should be noted that Tpetra_LAPACK_wrappers.h is an almost 
    exact duplicate of Epetra_LAPACK_wrappers.h. 
    More info on how/why to be added 
    at a later date, when Paul gets around to it, and once he has a better idea 
    of what's going on here.
*/
/* Define fcd (Fortran Tpetra_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)


#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define Tpetra_fcd fcd 


#elif defined(INTEL_CXML)

#define PREFIX __stdcall 
#define Tpetra_fcd char *, unsigned int 


#elif defined(INTEL_MKL)

#define PREFIX
#define Tpetra_fcd char *

#endif 

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_FUNC is defined undefine it because we want to redefine */

#ifdef F77_FUNC
#undef F77_FUNC
#endif


#define F77_FUNC(lcase,UCASE) PREFIX UCASE

#else /* Define Tpetra_fcd for all other machines */

#define PREFIX
#define Tpetra_fcd char * 

/* In the future use autoconf's definition of F77_FUNC */ 
#ifdef F77_FUNC
#undef F77_FUNC
#endif

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(lcase,UCASE) lcase ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(lcase,UCASE) lcase ## __

#endif


#define DGETRF_F77  F77_FUNC(dgetrf,DGETRF)
#define DGETRS_F77  F77_FUNC(dgetrs,DGETRS)
#define DGETRI_F77  F77_FUNC(dgetri,DGETRI)
#define DGERFS_F77  F77_FUNC(dgerfs,DGERFS)
#define DGECON_F77  F77_FUNC(dgecon,DGECON)
#define DGESVX_F77  F77_FUNC(dgesvx,DGESVX)
#define DGESV_F77   F77_FUNC(dgesv,DGESV)
#define DGEEQU_F77  F77_FUNC(dgeequ,DGEEQU)
#define DPOTRF_F77  F77_FUNC(dpotrf,DPOTRF)
#define DPOTRS_F77  F77_FUNC(dpotrs,DPOTRS)
#define DPOTRI_F77  F77_FUNC(dpotri,DPOTRI)
#define DPOCON_F77  F77_FUNC(dpocon,DPOCON)
#define DPOSV_F77   F77_FUNC(dposv,DPOSV)
#define DPOEQU_F77  F77_FUNC(dpoequ,DPOEQU)
#define DPORFS_F77  F77_FUNC(dporfs,DPORFS)
#define DPOSVX_F77  F77_FUNC(dposvx,DPOSVX)
#define DLAMCH_F77  F77_FUNC(dlamch,DLAMCH)
#define DGELS_F77   F77_FUNC(dgels,DGELS)
#define DGEEV_F77   F77_FUNC(dgeev,DGEEV)
#define DGEHRD_F77  F77_FUNC(dgehrd,DGEHRD)
#define DHSEQR_F77  F77_FUNC(dhseqr,DHSEQR)
#define DORGHR_F77  F77_FUNC(dorghr,DORGHR)
#define DORMHR_F77  F77_FUNC(dormhr,DORMHR)
#define DTREVC_F77  F77_FUNC(dtrevc,DTREVC)
#define DTREXC_F77  F77_FUNC(dtrexc,DTREXC)

#define SGETRF_F77  F77_FUNC(sgetrf,SGETRF)
#define SGETRS_F77  F77_FUNC(sgetrs,SGETRS)
#define SGETRI_F77  F77_FUNC(sgetri,SGETRI)
#define SGERFS_F77  F77_FUNC(sgerfs,SGERFS)
#define SGECON_F77  F77_FUNC(sgecon,SGECON)
#define SGESVX_F77  F77_FUNC(sgesvx,SGESVX)
#define SGESV_F77   F77_FUNC(sgesv,SGESV)
#define SGEEQU_F77  F77_FUNC(sgeequ,SGEEQU)
#define SPOTRF_F77  F77_FUNC(spotrf,SPOTRF)
#define SPOTRS_F77  F77_FUNC(spotrs,SPOTRS)
#define SPOTRI_F77  F77_FUNC(spotri,SPOTRI)
#define SPOCON_F77  F77_FUNC(spocon,SPOCON)
#define SPOSV_F77   F77_FUNC(sposv,SPOSV)
#define SPOEQU_F77  F77_FUNC(spoequ,SPOEQU)
#define SPORFS_F77  F77_FUNC(sporfs,SPORFS)
#define SPOSVX_F77  F77_FUNC(sposvx,SPOSVX)
#define SGELS_F77   F77_FUNC(sgels,SGELS)
#define SGEEV_F77   F77_FUNC(sgeev,SGEEV)
#define SGEHRD_F77  F77_FUNC(sgehrd,SGEHRD)
#define SHSEQR_F77  F77_FUNC(shseqr,SHSEQR)
#define SORGHR_F77  F77_FUNC(sorghr,SORGHR)
#define SORMHR_F77  F77_FUNC(sormhr,SORMHR)
#define STREVC_F77  F77_FUNC(strevc,STREVC)
#define STREXC_F77  F77_FUNC(strexc,STREXC)
#define SLAMCH_F77  F77_FUNC(slamch,SLAMCH)

#ifdef __cplusplus
extern "C" {
#endif

// Double precision LAPACK linear solvers
void DGETRF_F77(int* m, int* n, double* a, int* lda, int* ipiv, int* info); 
void DGETRS_F77(Tpetra_fcd, int* n, int* nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void DGETRI_F77(int* n, double* a, int* lda, int*ipiv, double * work , int* lwork, int* info);
void DGECON_F77(Tpetra_fcd norm, int* n, double* a, int* lda, 
                       double *anorm, double * rcond, double * work,
                       int * iwork, int* info); 
void DGESV_F77(int * n, int * nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void DGEEQU_F77(int* m, int* n, double* a, int* lda, double * r, double * c, 
			double * rowcnd, double * colcnd,
                       double * amax, int* info); 
void DGERFS_F77(Tpetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, int*ipiv, 
                       double * b, int * ldb, double * x, int * ldx, 
                       double * ferr, double * berr, double * work, 
                       int * iwork, int * info);
void DGESVX_F77(Tpetra_fcd, Tpetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, int*ipiv, Tpetra_fcd, 
                       double * r, double *c, double * b, int * ldb, double * x, int * ldx, 
                       double * rcond, double * ferr, double * berr, double * work, 
                       int * iwork, int * info);

void DPOTRF_F77(Tpetra_fcd, int* n, double* a, int* lda, int* info); 
void DPOTRS_F77(Tpetra_fcd, int * n, int * nrhs, double* a,
                       int* lda, double*x , int* ldx, int* info);
void DPOTRI_F77(Tpetra_fcd, int* n, double* a, int* lda, int* info); 
void DPOCON_F77(Tpetra_fcd, int* n, double* a, int* lda, 
                       double * anorm, double * rcond, double * work,
                       int * iwork, int* info); 
void DPOSV_F77(Tpetra_fcd, int * n, int * nrhs, double* a,
                       int* lda, double*x , int* ldx, int* info);
void DPOEQU_F77(int* n, double* a, int* lda, double * s, double * scond,
                       double * amax, int* info); 

void DPORFS_F77(Tpetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, 
                       double * b, int * ldb, double * x, int * ldx, 
                       double * ferr, double * berr, double * work, 
                       int * iwork, int * info);

void DPOSVX_F77(Tpetra_fcd, Tpetra_fcd, int * n, int * nrhs, double * a, 
                       int * lda, double * af, int * ldaf, Tpetra_fcd, 
                       double * s, double * b, int * ldb, double * x, int * ldx, 
                       double * rcond, double * ferr, double * berr, double * work, 
                       int * iwork, int * info);

// Single precision LAPACK linear solvers
void SGETRF_F77(int* m, int* n, float* a, int* lda, int* ipiv, int* info); 
void SGETRS_F77(Tpetra_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
void SGETRI_F77(int* n, float* a, int* lda, int*ipiv, float * work , int* lwork, int* info);
void SGECON_F77(Tpetra_fcd norm, int* n, float* a, int* lda, 
                       float * anorm, float * rcond, float * work,
                       int * iwork, int* info); 
void SGESV_F77(int * n, int * nrhs, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
void SGEEQU_F77(int* m, int* n, float* a, int* lda, float * r, float * c, 
			float * rowcnd, float * colcnd,
			float * amax, int* info); 
void SGERFS_F77(Tpetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, int*ipiv, 
                       float * b, int * ldb, float * x, int * ldx, 
                       float * ferr, float * berr, float * work, 
                       int * iwork, int * info);
void SGESVX_F77(Tpetra_fcd, Tpetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, int*ipiv, Tpetra_fcd, 
                       float * r, float *c, float * b, int * ldb, float * x, int * ldx, 
                       float * rcond, float * ferr, float * berr, float * work, 
                       int * iwork, int * info);

void SPOTRF_F77(Tpetra_fcd, int* n, float* a, int* lda, int* info); 
void SPOTRS_F77(Tpetra_fcd, int * n, int * nrhs, float* a,
                       int* lda, float*x , int* ldx, int* info);
void SPOTRI_F77(Tpetra_fcd, int* n, float* a, int* lda, int* info); 
void SPOCON_F77(Tpetra_fcd, int* n, float* a, int* lda, 
                       float * anorm, float * rcond, float * work,
                       int * iwork, int* info); 
void SPOSV_F77(Tpetra_fcd, int * n, int * nrhs, float* a,
                       int* lda, float*x , int* ldx, int* info);
void SPOEQU_F77(int* n, float* a, int* lda, float * s, float * scond,
                       float * amax, int* info); 

void SPORFS_F77(Tpetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, 
                       float * b, int * ldb, float * x, int * ldx, 
                       float * ferr, float * berr, float * work, 
                       int * iwork, int * info);

void SPOSVX_F77(Tpetra_fcd, Tpetra_fcd, int * n, int * nrhs, float * a, 
                       int * lda, float * af, int * ldaf, Tpetra_fcd, 
                       float * s, float * b, int * ldb, float * x, int * ldx, 
                       float * rcond, float * ferr, float * berr, float * work, 
                       int * iwork, int * info);

// Double precision LAPACK eigen solvers
void DGELS_F77(Tpetra_fcd ch, int*, int*, int*,
                       double*, int*, double*, int*, double*, int*, int*);

void DGEEV_F77(Tpetra_fcd, Tpetra_fcd, int*, double*, int*,
                      double*, double*, double*, int*, double*, int*,
                      double*, int*, int*);
           
void  DGEHRD_F77(int * n, int * ilo, int * ihi, double * A,
                        int * lda, double * tau, double * work, int * lwork,
                        int * info);

void  DHSEQR_F77(Tpetra_fcd job, Tpetra_fcd, int * n, int * ilo, int * ihi,
                        double * h, int * ldh, double * wr, double * wi, double * z,
                        int * ldz, double * work, int * lwork, int * info);

void  DORGHR_F77(int * n, int * ilo, int * ihi, double * a, int * lda, double * tau,
                        double * work, int * lwork, int * info);
                        
void  DORMHR_F77(Tpetra_fcd, Tpetra_fcd, int * m, int * n, int * ilo,
                        int * ihi, double * a, int * lda, double * tau, double * c,
                        int * ldc, double * work, int * lwork, int * info);

void  DTREVC_F77(Tpetra_fcd, Tpetra_fcd, int * select, int * n, double * t,
                        int * ldt, double *vl, int * ldvl, double * vr, int * ldvr,
                        int * mm, int * m, double * work, int * info); 

void  DTREXC_F77(Tpetra_fcd, int * n, double * t, int * ldt, double * q,
                        int * ldq, int * ifst, int * ilst, double * work, int * info);

double  DLAMCH_F77(Tpetra_fcd);


// Single precision LAPACK eigen solvers
void SGELS_F77(Tpetra_fcd, int*, int*, int*,
                       float*, int*, float*, int*, float*, int*, int*);

void SGEEV_F77(Tpetra_fcd, Tpetra_fcd, int*, float*, int*,
                      float*, float*, float*, int*, float*, int*,
                      float*, int*, int*);

void  SGEHRD_F77(int * n, int * ilo, int * ihi, float * A,
                        int * lda, float * tau, float * work, int * lwork,
                        int * info);

void  SHSEQR_F77(Tpetra_fcd job, Tpetra_fcd, int * n, int * ilo, int * ihi,
                        float * h, int * ldh, float * wr, float * wi, float * z,
                        int * ldz, float * work, int * lwork, int * info);

void  SORGHR_F77(int * n, int * ilo, int * ihi, float * a, int * lda, float * tau,
                        float * work, int * lwork, int * info);
                        
void  SORMHR_F77(Tpetra_fcd, Tpetra_fcd, int * m, int * n, int * ilo,
                        int * ihi, float * a, int * lda, float * tau, float * c,
                        int * ldc, float * work, int * lwork, int * info);

void  STREVC_F77(Tpetra_fcd, Tpetra_fcd, int * select, int * n, float * t,
                        int * ldt, float *vl, int * ldvl, float * vr, int * ldvr,
                        int * mm, int * m, float * work, int * info); 

void  STREXC_F77(Tpetra_fcd, int * n, float * t, int * ldt, float * q,
                        int * ldq, int * ifst, int * ilst, float * work, int * info);

float  SLAMCH_F77(Tpetra_fcd);

#ifdef __cplusplus
}
#endif

#endif // end of TPETRA_LAPACK_WRAPPERS_H_
