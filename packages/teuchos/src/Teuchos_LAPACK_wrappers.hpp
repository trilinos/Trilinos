// Kris
// 06.11.03 -- Format cleanup
// 06.17.03 -- Added LAPY2 and GEES by request
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_LAPACK_WRAPPERS_HPP_
#define _TEUCHOS_LAPACK_WRAPPERS_HPP_

#include "Teuchos_ConfigDefs.hpp"

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
#define Teuchos_fcd fcd 

#define DGETRF_F77  F77_FUNC(sgetrf,SGETRF)
#define DGETRS_F77  F77_FUNC(sgetrs,SGETRS)
#define DGETRI_F77  F77_FUNC(sgetri,SGETRI)
#define DGERFS_F77  F77_FUNC(sgerfs,SGERFS)
#define DGECON_F77  F77_FUNC(sgecon,SGECON)
#define DGESVX_F77  F77_FUNC(sgesvx,SGESVX)
#define DGESV_F77   F77_FUNC(sgesv,SGESV)
#define DGEEQU_F77  F77_FUNC(sgeequ,SGEEQU)
#define DPOTRF_F77  F77_FUNC(spotrf,SPOTRF)
#define DPOTRS_F77  F77_FUNC(spotrs,SPOTRS)
#define DPOTRI_F77  F77_FUNC(dpotri,SPOTRI)
#define DPOCON_F77  F77_FUNC(spocon,SPOCON)
#define DPOSV_F77   F77_FUNC(sposv,SPOSV)
#define DPOEQU_F77  F77_FUNC(spoequ,SPOEQU)
#define DPORFS_F77  F77_FUNC(sporfs,SPORFS)
#define DPOSVX_F77  F77_FUNC(sposvx,SPOSVX)
#define DLAMCH_F77  F77_FUNC(slamch,SLAMCH)
#define DGELS_F77   F77_FUNC(sgels,SGELS)
#define DGEEV_F77   F77_FUNC(sgeev,SGEEV)
#define DGEHRS_F77  F77_FUNC(sgehrd,SGEHRD)
#define DHSEQR_F77  F77_FUNC(shseqr,SHSEQR)
#define DORGHR_F77  F77_FUNC(sorghr,SORGHR)
#define DORMHR_F77  F77_FUNC(sormhr,SORMHR)
#define DTREVC_F77  F77_FUNC(strevc,STREVC)
#define DTREXC_F77  F77_FUNC(strexc,STREXC)
#define DGEES_F77   F77_FUNC(sgees,SGEES)
#define DLAPY2_F77  F77_FUNC(slapy2,SLAPY2)

#elif defined(INTEL_CXML)

#define PREFIX __stdcall 
#define Teuchos_fcd char *, unsigned int 

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
#define DGEES_F77   F77_FUNC(dgees,DGEES)
#define DLAPY2_F77  F77_FUNC(dlapy2,DLAPY2)

#elif defined(INTEL_MKL)

#define PREFIX
#define Teuchos_fcd char *

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
#define DGEES_F77   F77_FUNC(dgees,DGEES)
#define DLAPY2_F77  F77_FUNC(dlapy2,DLAPY2)

#endif 

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_FUNC is defined undefine it because we want to redefine */

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#define F77_FUNC(lcase,UCASE) PREFIX UCASE

#else /* Define Teuchos_fcd for all other machines */

#define PREFIX
#define Teuchos_fcd char * 

#ifndef HAVE_CONFIG_H

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_FUNC(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_FUNC(lcase,UCASE) lcase ## _
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */

#endif /* HAVE_CONFIG_H */

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
#define DGEES_F77   F77_FUNC(dgees,DGEES)
#define DLAPY2_F77  F77_FUNC(dlapy2,DLAPY2)

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
#define SGEES_F77   F77_FUNC(sgees,SGEES)
#define SLAPY2_F77  F77_FUNC(slapy2,SLAPY2)


#ifdef __cplusplus
extern "C" {
#endif

// Double precision LAPACK linear solvers
void PREFIX DGETRF_F77(int* m, int* n, double* a, int* lda, int* ipiv, int* info); 
void PREFIX DGETRS_F77(Teuchos_fcd, int* n, int* nrhs, double* a, int* lda, int*ipiv, double*x , int* ldx, int* info);
void PREFIX DGETRI_F77(int* n, double* a, int* lda, int*ipiv, double* work , int* lwork, int* info);
void PREFIX DGECON_F77(Teuchos_fcd norm, int* n, double* a, int* lda, double* anorm, double* rcond, double* work, int* iwork, int* info); 
void PREFIX DGESV_F77(int* n, int* nrhs, double* a, int* lda, int*ipiv, double*x , int* ldx, int* info);
void PREFIX DGEEQU_F77(int* m, int* n, double* a, int* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, int* info); 
void PREFIX DGERFS_F77(Teuchos_fcd, int* n, int* nrhs, double* a, int* lda, double* af, int* ldaf, int*ipiv, double* b, int* ldb, double* x, int* ldx, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DGESVX_F77(Teuchos_fcd, Teuchos_fcd, int* n, int* nrhs, double* a, int* lda, double* af, int* ldaf, int*ipiv, Teuchos_fcd, double* r,
double* c, double* b, int* ldb, double* x, int* ldx, double* rcond, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DPOTRF_F77(Teuchos_fcd, int* n, double* a, int* lda, int* info); 
void PREFIX DPOTRS_F77(Teuchos_fcd, int* n, int* nrhs, double* a, int* lda, double*x , int* ldx, int* info);
void PREFIX DPOTRI_F77(Teuchos_fcd, int* n, double* a, int* lda, int* info); 
void PREFIX DPOCON_F77(Teuchos_fcd, int* n, double* a, int* lda, double* anorm, double* rcond, double* work, int* iwork, int* info); 
void PREFIX DPOSV_F77(Teuchos_fcd, int* n, int* nrhs, double* a, int* lda, double*x , int* ldx, int* info);
void PREFIX DPOEQU_F77(int* n, double* a, int* lda, double* s, double* scond, double* amax, int* info); 
void PREFIX DPORFS_F77(Teuchos_fcd, int* n, int* nrhs, double* a, int* lda, double* af, int* ldaf, double* b, int* ldb, double* x, int* ldx, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DPOSVX_F77(Teuchos_fcd, Teuchos_fcd, int* n, int* nrhs, double* a, int* lda, double* af, int* ldaf, Teuchos_fcd, double* s, double* b, int* ldb, double* x, int* ldx, double* rcond, double* ferr, double* berr, double* work, int* iwork, int* info);

// Single precision LAPACK linear solvers
void PREFIX SGETRF_F77(int* m, int* n, float* a, int* lda, int* ipiv, int* info); 
void PREFIX SGETRS_F77(Teuchos_fcd, int* m, int* n, float* a, int* lda, int*ipiv, float*x , int* ldx, int* info);
void PREFIX SGETRI_F77(int* n, float* a, int* lda, int*ipiv, float* work , int* lwork, int* info);
void PREFIX SGECON_F77(Teuchos_fcd norm, int* n, float* a, int* lda, float* anorm, float* rcond, float* work, int* iwork, int* info); 
void PREFIX SGESV_F77(int* n, int* nrhs, float* a, int* lda, int*ipiv, float*x , int* ldx, int* info);
void PREFIX SGEEQU_F77(int* m, int* n, float* a, int* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, int* info); 
void PREFIX SGERFS_F77(Teuchos_fcd, int* n, int* nrhs, float* a, int* lda, float* af, int* ldaf, int*ipiv, float* b, int* ldb, float* x, int* ldx, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SGESVX_F77(Teuchos_fcd, Teuchos_fcd, int* n, int* nrhs, float* a, int* lda, float* af, int* ldaf, int*ipiv, Teuchos_fcd, float* r, float* c, float* b, int* ldb, float* x, int* ldx, float* rcond, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SPOTRF_F77(Teuchos_fcd, int* n, float* a, int* lda, int* info); 
void PREFIX SPOTRS_F77(Teuchos_fcd, int* n, int* nrhs, float* a, int* lda, float*x , int* ldx, int* info);
void PREFIX SPOTRI_F77(Teuchos_fcd, int* n, float* a, int* lda, int* info); 
void PREFIX SPOCON_F77(Teuchos_fcd, int* n, float* a, int* lda, float* anorm, float* rcond, float* work, int* iwork, int* info); 
void PREFIX SPOSV_F77(Teuchos_fcd, int* n, int* nrhs, float* a, int* lda, float*x , int* ldx, int* info);
void PREFIX SPOEQU_F77(int* n, float* a, int* lda, float* s, float* scond, float* amax, int* info); 
void PREFIX SPORFS_F77(Teuchos_fcd, int* n, int* nrhs, float* a, int* lda, float* af, int* ldaf, float* b, int* ldb, float* x, int* ldx, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SPOSVX_F77(Teuchos_fcd, Teuchos_fcd, int* n, int* nrhs, float* a, int* lda, float* af, int* ldaf, Teuchos_fcd, float* s, float* b, int* ldb, float* x, int* ldx, float* rcond, float* ferr, float* berr, float* work, int* iwork, int* info);

// Double precision LAPACK eigen solvers
void PREFIX DGELS_F77(Teuchos_fcd ch, int*, int*, int*, double*, int*, double*, int*, double*, int*, int*);
void PREFIX DGEEV_F77(Teuchos_fcd, Teuchos_fcd, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
void PREFIX DGEHRD_F77(int* n, int* ilo, int* ihi, double* A, int* lda, double* tau, double* work, int* lwork, int* info);
void PREFIX DHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, int* n, int* ilo, int* ihi, double* h, int* ldh, double* wr, double* wi, double* z, int* ldz, double* work, int* lwork, int* info);
void PREFIX DORGHR_F77(int* n, int* ilo, int* ihi, double* a, int* lda, double* tau, double* work, int* lwork, int* info);
void PREFIX DORMHR_F77(Teuchos_fcd, Teuchos_fcd, int* m, int* n, int* ilo, int* ihi, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* lwork, int* info);
void PREFIX DTREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, int* n, double* t, int* ldt, double* vl, int* ldvl, double* vr, int* ldvr, int* mm, int* m, double* work, int* info); 
void PREFIX DTREXC_F77(Teuchos_fcd, int* n, double* t, int* ldt, double* q, int* ldq, int* ifst, int* ilst, double* work, int* info);
double DLAMCH_F77(Teuchos_fcd);

// Single precision LAPACK eigen solvers
void PREFIX SGELS_F77(Teuchos_fcd, int*, int*, int*, float*, int*, float*, int*, float*, int*, int*);
void PREFIX SGEEV_F77(Teuchos_fcd, Teuchos_fcd, int*, float*, int*, float*, float*, float*, int*, float*, int*, float*, int*, int*);
void PREFIX SGEHRD_F77(int* n, int* ilo, int* ihi, float* A, int* lda, float* tau, float* work, int* lwork, int* info);
void PREFIX SHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, int* n, int* ilo, int* ihi, float* h, int* ldh, float* wr, float* wi, float* z, int* ldz, float* work, int* lwork, int* info);
void PREFIX SORGHR_F77(int* n, int* ilo, int* ihi, float* a, int* lda, float* tau, float* work, int* lwork, int* info);
void PREFIX SORMHR_F77(Teuchos_fcd, Teuchos_fcd, int* m, int* n, int* ilo, int* ihi, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info);
void PREFIX STREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, int* n, float* t, int* ldt, float* vl, int* ldvl, float* vr, int* ldvr, int* mm, int* m, float* work, int* info); 
void PREFIX STREXC_F77(Teuchos_fcd, int* n, float* t, int* ldt, float* q, int* ldq, int* ifst, int* ilst, float* work, int* info);
float SLAMCH_F77(Teuchos_fcd);

void PREFIX SGEES_F77(Teuchos_fcd, Teuchos_fcd, int*, int*, float*, int*, int*, float*, float*, float*, int*, float*, int*, int*, int*);
void PREFIX DGEES_F77(Teuchos_fcd, Teuchos_fcd, int*, int*, double*, int*, int*, double*, double*, double*, int*, double*, int*, int*, int*);
float PREFIX SLAPY2_F77(float*, float*);
double PREFIX DLAPY2_F77(double*, double*);

#ifdef __cplusplus
}

#endif

#endif // end of TEUCHOS_LAPACK_WRAPPERS_HPP_
