// Kris
// 06.11.03 -- Format cleanup
// 06.17.03 -- Added LAPY2 and GEES by request
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_LAPACK_WRAPPERS_HPP_
#define _TEUCHOS_LAPACK_WRAPPERS_HPP_

#include "Teuchos_ConfigDefs.hpp"

/*! \file Teuchos_LAPACK_wrappers.hpp
    \brief Document this!
*/
/* Define fcd (Fortran Teuchos_fcd descriptor) for non-standard situations */

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
#define DPOTRI_F77  F77_FUNC(spotri,SPOTRI)
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
#define Teuchos_fcd const char *, unsigned int 

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
#define Teuchos_fcd const char *

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
#define Teuchos_fcd const char * 

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
void PREFIX DGELS_F77(Teuchos_fcd ch, const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb, double* work, const int* lwork, int* info);
void PREFIX DGETRF_F77(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info); 
void PREFIX DGETRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda,const int* ipiv, double* x , const int* ldx, int* info);
void PREFIX DGETRI_F77(const int* n, double* a, const int* lda, const int* ipiv, double* work , const int* lwork, int* info);
void PREFIX DGECON_F77(Teuchos_fcd norm, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info); 
void PREFIX DGESV_F77(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* x , const int* ldx, int* info);
void PREFIX DGEEQU_F77(const int* m, const int* n, const double* a, const int* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, int* info); 
void PREFIX DGERFS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda, const double* af, const int* ldaf, const int* ipiv, const double* b, const int* ldb, double* x, const int* ldx, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DGESVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, double* af, const int* ldaf, int* ipiv, Teuchos_fcd, double* r,
double* c, double* b, const int* ldb, double* x, const int* ldx, double* rcond, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DPOTRF_F77(Teuchos_fcd, const int* n, double* a, const int* lda, int* info); 
void PREFIX DPOTRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda, double*x , const int* ldx, int* info);
void PREFIX DPOTRI_F77(Teuchos_fcd, const int* n, double* a, const int* lda, int* info); 
void PREFIX DPOCON_F77(Teuchos_fcd, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info); 
void PREFIX DPOSV_F77(Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, double*x , const int* ldx, int* info);
void PREFIX DPOEQU_F77(const int* n, const double* a, const int* lda, double* s, double* scond, double* amax, int* info); 
void PREFIX DPORFS_F77(Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, const double* af, const int* ldaf, const double* b, const int* ldb, double* x, const int* ldx, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DPOSVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, double* af, const int* ldaf, Teuchos_fcd, double* s, double* b, const int* ldb, double* x, const int* ldx, double* rcond, double* ferr, double* berr, double* work, int* iwork, int* info);

// Single precision LAPACK linear solvers
void PREFIX SGELS_F77(Teuchos_fcd ch, const int* m, const int* n, const int* nrhs, float* a, const int* lda, float* b, const int* ldb, float* work, const int* lwork, int* info);
void PREFIX SGETRF_F77(const int* m, const int* n, float* a, const int* lda, int* ipiv, int* info);void PREFIX SGETRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const float* a, const int* lda,const int* ipiv, float* x , const int* ldx, int* info);
void PREFIX SGETRI_F77(const int* n, float* a, const int* lda, const int* ipiv, float* work , const int* lwork, int* info);
void PREFIX SGECON_F77(Teuchos_fcd norm, const int* n, const float* a, const int* lda, const float* anorm, float* rcond, float* work, int* iwork, int* info); 
void PREFIX SGESV_F77(const int* n, const int* nrhs, float* a, const int* lda, int* ipiv, float* x , const int* ldx, int* info);
void PREFIX SGEEQU_F77(const int* m, const int* n, const float* a, const int* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, int* info); 
void PREFIX SGERFS_F77(Teuchos_fcd, const int* n, const int* nrhs, const float* a, const int* lda, const float* af, const int* ldaf, const int* ipiv, const float* b, const int* ldb, float* x, const int* ldx, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SGESVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, float* af, const int* ldaf, int* ipiv, Teuchos_fcd, float* r,
float* c, float* b, const int* ldb, float* x, const int* ldx, float* rcond, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SPOTRF_F77(Teuchos_fcd, const int* n, float* a, const int* lda, int* info); 
void PREFIX SPOTRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const float* a, const int* lda, float*x , const int* ldx, int* info);
void PREFIX SPOTRI_F77(Teuchos_fcd, const int* n, float* a, const int* lda, int* info); 
void PREFIX SPOCON_F77(Teuchos_fcd, const int* n, const float* a, const int* lda, const float* anorm, float* rcond, float* work, int* iwork, int* info); 
void PREFIX SPOSV_F77(Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, float*x , const int* ldx, int* info);
void PREFIX SPOEQU_F77(const int* n, const float* a, const int* lda, float* s, float* scond, float* amax, int* info); 
void PREFIX SPORFS_F77(Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, const float* af, const int* ldaf, const float* b, const int* ldb, float* x, const int* ldx, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SPOSVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, float* af, const int* ldaf, Teuchos_fcd, float* s, float* b, const int* ldb, float* x, const int* ldx, float* rcond, float* ferr, float* berr, float* work, int* iwork, int* info);

// Double precision LAPACK eigen solvers
void PREFIX DGEEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info);
void PREFIX DGEHRD_F77(const int* n, const int* ilo, const int* ihi, double* A, const int* lda, double* tau, double* work, const int* lwork, int* info);
void PREFIX DHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, const int* n, const int* ilo, const int* ihi, double* h, const int* ldh, double* wr, double* wi, double* z, const int* ldz, double* work, const int* lwork, int* info);
void PREFIX DGEES_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, double* a, const int* lda, int*sdim, double* wr, double* wi, double* vs, const int* ldvs, double* work, const int* lwork, int* bwork, int* info);
void PREFIX DORGHR_F77(const int* n, const int* ilo, const int* ihi, double* a, const int* lda, double* tau, double* work, int* lwork, int* info);
void PREFIX DORMHR_F77(Teuchos_fcd, Teuchos_fcd, const int* m, const int* n, const int* ilo, const int* ihi, const double* a, const int* lda, const double* tau, double* c, const int* ldc, double* work, int* lwork, int* info);
void PREFIX DTREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, const double* t, const int* ldt, double* vl, const int* ldvl, double* vr, const int* ldvr, const int* mm, int* m, double* work, int* info); 
void PREFIX DTREXC_F77(Teuchos_fcd, const int* n, double* t, const int* ldt, double* q, const int* ldq, int* ifst, int* ilst, double* work, int* info);

// Single precision LAPACK eigen solvers
void PREFIX SGEEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, float* a, const int* lda, float* wr, float* wi, float* vl, const int* ldvl, float* vr, const int* ldvr, float* work, const int* lwork, int* info);
void PREFIX SGEHRD_F77(const int* n, const int* ilo, const int* ihi, float* A, const int* lda, float* tau, float* work, const int* lwork, int* info);
void PREFIX SHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, const int* n, const int* ilo, const int* ihi, float* h, const int* ldh, float* wr, float* wi, float* z, const int* ldz, float* work, const int* lwork, int* info);
void PREFIX SGEES_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, float* a, const int* lda, int* sdim, float* wr, float* wi, float* vs, const int* ldvs, float* work, const int* lwork, int* bwork, int* info);
void PREFIX SORGHR_F77(const int* n, const int* ilo, const int* ihi, float* a, const int* lda, float* tau, float* work, int* lwork, int* info);
void PREFIX SORMHR_F77(Teuchos_fcd, Teuchos_fcd, const int* m, const int* n, const int* ilo, const int* ihi, const float* a, const int* lda, const float* tau, float* c, const int* ldc, float* work, int* lwork, int* info);
void PREFIX STREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, const float* t, const int* ldt, float* vl, const int* ldvl, float* vr, const int* ldvr, const int* mm, int* m, float* work, int* info); 
void PREFIX STREXC_F77(Teuchos_fcd, const int* n, float* t, const int* ldt, float* q, const int* ldq, int* ifst, int* ilst, float* work, int* info);

float SLAMCH_F77(Teuchos_fcd);
double DLAMCH_F77(Teuchos_fcd);

float PREFIX SLAPY2_F77(const float* x, const float* y);
double PREFIX DLAPY2_F77(const double* x, const double* y);

#ifdef __cplusplus
}

#endif

#endif // end of TEUCHOS_LAPACK_WRAPPERS_HPP_
