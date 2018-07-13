/*
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
*/

#ifndef EPETRA_LAPACK_WRAPPERS_H
#define EPETRA_LAPACK_WRAPPERS_H

#include "Epetra_ConfigDefs.h"
#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)

#ifdef CRAY_T3X

#include "fortran.h"
#define Epetra_fcd fcd
#define PREFIX

/* CRAY Single precision is used like everyone else's double precision */
#define DGECON_F77  SGECON
#define DGEEQU_F77  SGEEQU
#define DGEEV_F77   SGEEV
#define DGEEVX_F77  SGEEVX
#define DGEHRD_F77  SGEHRD
#define DGELS_F77   SGELS
#define DGELSS_F77  SGELSS
#define DGEQPF_F77  SGEQPF
#define DGERFS_F77  SGERFS
#define DGESDD_F77  SGESDD
#define DGESVD_F77  SGESVD
#define DGESV_F77   SGESV
#define DGESVX_F77  SGESVX
#define DGETRF_F77  SGETRF
#define DGEQRF_F77  SGEQRF
#define DGETRI_F77  SGETRI
#define DGETRS_F77  SGETRS
#define DGGEV_F77   SGGEV
#define DGGLSE_F77  SGGLSE
#define DGGSVD_F77  SGGSVD
#define DHSEQR_F77  SHSEQR
#define DLAIC1_F77  SLAIC1
#define DLAMCH_F77  SLAMCH
#define DLARFT_F77  SLARFT
#define DLASWP_F77  SLASWP
#define DORGQR_F77  SORGQR
#define DORGHR_F77  SORGHR
#define DORMHR_F77  SORMHR
#define DPOCON_F77  SPOCON
#define DPOEQU_F77  SPOEQU
#define DPORFS_F77  SPORFS
#define DPOSV_F77   SPOSV
#define DPOSVX_F77  SPOSVX
#define DPOTRF_F77  SPOTRF
#define DPOTRI_F77  SPOTRI
#define DPOTRS_F77  SPOTRS
#define DSPEV_F77   SSPEV
#define DSPGV_F77   SSPGV
#define DSTEV_F77   SSTEV
#define DSYEVD_F77  SSYEVD
#define DSYEV_F77   SSYEV
#define DSYEVR_F77  SSYEVR
#define DSYEVX_F77  SSYEVX
#define DSYGV_F77   SSYGV
#define DSYGVX_F77  SSYGVX
#define DTREVC_F77  STREVC
#define DTREXC_F77  STREXC
#define DTRTRS_F77  STRTRS
/* Done with T3X double precision */
#endif

#if defined(INTEL_CXML)

#define Epetra_fcd const char *, const unsigned int
#define PREFIX __stdcall

#endif

#if defined(INTEL_MKL)

#define Epetra_fcd const char *
#define PREFIX

#endif

/* The remainder of this block is for T3X, CXML and MKL */

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif

#define F77_BLAS_MANGLE(lcase,UCASE) UCASE

#else
/* Not defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL) */

#define Epetra_fcd const char *
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
#endif /* defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL) */

#ifndef CRAY_T3X

#define DGECON_F77  F77_BLAS_MANGLE(dgecon,DGECON)
#define DGEEQU_F77  F77_BLAS_MANGLE(dgeequ,DGEEQU)
#define DGEEV_F77   F77_BLAS_MANGLE(dgeev,DGEEV)
#define DGEEVX_F77  F77_BLAS_MANGLE(dgeevx,DGEEVX)
#define DGEHRD_F77  F77_BLAS_MANGLE(dgehrd,DGEHRD)
#define DGELS_F77   F77_BLAS_MANGLE(dgels,DGELS)
#define DGELSS_F77  F77_BLAS_MANGLE(dgelss,DGELSS)
#define DGEQPF_F77  F77_BLAS_MANGLE(dgeqpf,DGEQPF)
#define DGERFS_F77  F77_BLAS_MANGLE(dgerfs,DGERFS)
#define DGESDD_F77  F77_BLAS_MANGLE(dgesdd,DGESDD)
#define DGESVD_F77  F77_BLAS_MANGLE(dgesvd,DGESVD)
#define DGESV_F77   F77_BLAS_MANGLE(dgesv,DGESV)
#define DGESVX_F77  F77_BLAS_MANGLE(dgesvx,DGESVX)
#define DGETRF_F77  F77_BLAS_MANGLE(dgetrf,DGETRF)
#define DGEQRF_F77  F77_BLAS_MANGLE(dgeqrf,DGEQRF)
#define DGETRI_F77  F77_BLAS_MANGLE(dgetri,DGETRI)
#define DGETRS_F77  F77_BLAS_MANGLE(dgetrs,DGETRS)
#define DGGEV_F77   F77_BLAS_MANGLE(dggev,DGGEV)
#define DGGLSE_F77  F77_BLAS_MANGLE(dgglse,DGGLSE)

#ifdef HAVE_EPETRA_LAPACK_GSSVD3
#define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd3,DGGSVD)
#else
#define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd,DGGSVD)
#endif

#define DHSEQR_F77  F77_BLAS_MANGLE(dhseqr,DHSEQR)
#define DLAIC1_F77  F77_BLAS_MANGLE(dlaic1,DLAIC1)
#define DLAMCH_F77  F77_BLAS_MANGLE(dlamch,DLAMCH)
#define DLARFT_F77  F77_BLAS_MANGLE(dlarft,DLARFT)
#define DLASWP_F77  F77_BLAS_MANGLE(dlaswp,DLASWP)
#define DORGQR_F77  F77_BLAS_MANGLE(dorgqr,DORGQR)
#define DORGHR_F77  F77_BLAS_MANGLE(dorghr,DORGHR)
#define DORMHR_F77  F77_BLAS_MANGLE(dormhr,DORMHR)
#define DPOCON_F77  F77_BLAS_MANGLE(dpocon,DPOCON)
#define DPOEQU_F77  F77_BLAS_MANGLE(dpoequ,DPOEQU)
#define DPORFS_F77  F77_BLAS_MANGLE(dporfs,DPORFS)
#define DPOSV_F77   F77_BLAS_MANGLE(dposv,DPOSV)
#define DPOSVX_F77  F77_BLAS_MANGLE(dposvx,DPOSVX)
#define DPOTRF_F77  F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define DPOTRI_F77  F77_BLAS_MANGLE(dpotri,DPOTRI)
#define DPOTRS_F77  F77_BLAS_MANGLE(dpotrs,DPOTRS)
#define DSPEV_F77   F77_BLAS_MANGLE(dspev,DSPEV)
#define DSPGV_F77   F77_BLAS_MANGLE(dspgv,DSPGV)
#define DSTEV_F77   F77_BLAS_MANGLE(dstev,DSTEV)
#define DSYEVD_F77  F77_BLAS_MANGLE(dsyevd,DSYEVD)
#define DSYEV_F77   F77_BLAS_MANGLE(dsyev,DSYEV)
#define DSYEVR_F77  F77_BLAS_MANGLE(dsyevr,DSYEVR)
#define DSYEVX_F77  F77_BLAS_MANGLE(dsyevx,DSYEVX)
#define DSYGV_F77   F77_BLAS_MANGLE(dsygv,DSYGV)
#define DSYGVX_F77  F77_BLAS_MANGLE(dsygvx,DSYGVX)
#define DTREVC_F77  F77_BLAS_MANGLE(dtrevc,DTREVC)
#define DTREXC_F77  F77_BLAS_MANGLE(dtrexc,DTREXC)
#define DTRTRS_F77  F77_BLAS_MANGLE(dtrtrs,DTRTRS)

/* End of defines for double precision when not on a T3X */

#endif

/* The following defines are good for all platforms */

#define SGECON_F77  F77_BLAS_MANGLE(sgecon,SGECON)
#define SGEEQU_F77  F77_BLAS_MANGLE(sgeequ,SGEEQU)
#define SGEEV_F77   F77_BLAS_MANGLE(sgeev,SGEEV)
#define SGEEVX_F77  F77_BLAS_MANGLE(sgeevx,SGEEVX)
#define SGEHRD_F77  F77_BLAS_MANGLE(sgehrd,SGEHRD)
#define SGELS_F77   F77_BLAS_MANGLE(sgels,SGELS)
#define SGELSS_F77  F77_BLAS_MANGLE(sgelss,SGELSS)
#define SGEQPF_F77  F77_BLAS_MANGLE(sgeqpf,SGEQPF)
#define SGERFS_F77  F77_BLAS_MANGLE(sgerfs,SGERFS)
#define SGESDD_F77  F77_BLAS_MANGLE(sgesdd,SGESDD)
#define SGESVD_F77  F77_BLAS_MANGLE(sgesvd,SGESVD)
#define SGESV_F77   F77_BLAS_MANGLE(sgesv,SGESV)
#define SGESVX_F77  F77_BLAS_MANGLE(sgesvx,SGESVX)
#define SGETRF_F77  F77_BLAS_MANGLE(sgetrf,SGETRF)
#define SGEQRF_F77  F77_BLAS_MANGLE(sgeqrf,SGEQRF)
#define SGETRI_F77  F77_BLAS_MANGLE(sgetri,SGETRI)
#define SGETRS_F77  F77_BLAS_MANGLE(sgetrs,SGETRS)
#define SGGEV_F77   F77_BLAS_MANGLE(sggev,SGGEV)
#define SGGLSE_F77  F77_BLAS_MANGLE(sgglse,SGGLSE)

#ifdef HAVE_EPETRA_LAPACK_GSSVD3
#define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd3,SGGSVD)
#else
#define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd,SGGSVD)
#endif

#define SHSEQR_F77  F77_BLAS_MANGLE(shseqr,SHSEQR)
#define SLAMCH_F77  F77_BLAS_MANGLE(slamch,SLAMCH)
#define SLARFT_F77  F77_BLAS_MANGLE(slarft,SLARFT)
#define SORGQR_F77  F77_BLAS_MANGLE(sorgqr,SORGQR)
#define SORGHR_F77  F77_BLAS_MANGLE(sorghr,SORGHR)
#define SORMHR_F77  F77_BLAS_MANGLE(sormhr,SORMHR)
#define SPOCON_F77  F77_BLAS_MANGLE(spocon,SPOCON)
#define SPOEQU_F77  F77_BLAS_MANGLE(spoequ,SPOEQU)
#define SPORFS_F77  F77_BLAS_MANGLE(sporfs,SPORFS)
#define SPOSV_F77   F77_BLAS_MANGLE(sposv,SPOSV)
#define SPOSVX_F77  F77_BLAS_MANGLE(sposvx,SPOSVX)
#define SPOTRF_F77  F77_BLAS_MANGLE(spotrf,SPOTRF)
#define SPOTRI_F77  F77_BLAS_MANGLE(spotri,SPOTRI)
#define SPOTRS_F77  F77_BLAS_MANGLE(spotrs,SPOTRS)
#define SSPEV_F77   F77_BLAS_MANGLE(sspev,SSPEV)
#define SSPGV_F77   F77_BLAS_MANGLE(sspgv,SSPGV)
#define SSTEV_F77   F77_BLAS_MANGLE(sstev,SSTEV)
#define SSYEVD_F77  F77_BLAS_MANGLE(ssyevd,SSYEVD)
#define SSYEV_F77   F77_BLAS_MANGLE(ssyev,SSYEV)
#define SSYEVR_F77  F77_BLAS_MANGLE(ssyevr,SSYEVR)
#define SSYEVX_F77  F77_BLAS_MANGLE(ssyevx,SSYEVX)
#define SSYGV_F77   F77_BLAS_MANGLE(ssygv,SSYGV)
#define SSYGVX_F77  F77_BLAS_MANGLE(ssygvx,SSYGVX)
#define STREVC_F77  F77_BLAS_MANGLE(strevc,STREVC)
#define STREXC_F77  F77_BLAS_MANGLE(strexc,STREXC)
#define STRTRS_F77  F77_BLAS_MANGLE(strtrs,STRTRS)

#ifdef __cplusplus
extern "C" {
#endif


  void PREFIX DGECON_F77(Epetra_fcd norm, const int* n, const double* a, const int* lda, const double *anorm, double * rcond,
			 double * work, int * iwork, int* info);
  void PREFIX DGEEQU_F77(const int* m, const int* n, const double* a, const int* lda, double * r, double * c, double * rowcnd,
			 double * colcnd, double * amax, int* info);
  void PREFIX DGEEV_F77(Epetra_fcd, Epetra_fcd, const int* n, double* a, const int* lda, double* wr, double* wi,
			double* vl, const int* ldvl,
			double* vr, const int* ldvr, double* work, const int* lwork, int* info);
  void PREFIX DGEEVX_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, double * a, const int * lda,
			 double * wr, double * wi, double * vl, const int * ldvl, double * vr, const int * ldvr,
			 int * ilo, int * ihi, double * scale, double * abnrm, double * rconde, double * rcondv,
			 double * work, const int * lwork, int * iwork, int * info);
  void PREFIX DGEHRD_F77(const int * n, const int * ilo, const int * ihi, double * A, const int * lda, double * tau, double * work,
			 const int * lwork, int * info);
  void PREFIX DGELS_F77(Epetra_fcd ch, const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb,
			double* work, const int* lwork, int* info);
  void PREFIX DGELSS_F77(const int * m, const int * n, const int * nrhs, double * a, const int * lda, double * b, const int * ldb,
			 double * s, const double * rcond, int * rank, double * work, const int * lwork, int * info);
  void PREFIX DGEQPF_F77(const int * m, const int * n, double * a, const int * lda, int * jpvt, double * tau, double * work, int * info);
  void PREFIX DGERFS_F77(Epetra_fcd, const int * n, const int * nrhs, const double * a, const int * lda, const double * af, const int * ldaf,
			 const int*ipiv, const double * b, const int * ldb, double * x, const int * ldx, double * ferr, double * berr,
			 double * work, int * iwork, int * info);

  void PREFIX DGESDD_F77(Epetra_fcd, const int * m, const int * n, double * a, const int * lda, double * s, double * u,
			 const int * ldu, double * vt, const int * ldvt, double * work, const int * lwork, int * iwork, int * info);

  void PREFIX DGESVD_F77(Epetra_fcd, Epetra_fcd, const int* m, const int* n, double* a, const int* lda, double* s, double* u,
			 const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info);
  void PREFIX DGESV_F77(const int * n, const int * nrhs, double* a, const int* lda, int*ipiv, double*x , const int* ldx, int* info);
  void PREFIX DGESVX_F77(Epetra_fcd, Epetra_fcd, const int * n, const int * nrhs, double * a, const int * lda, double * af,
			 const int * ldaf, int*ipiv, Epetra_fcd, double * r, double *c, double * b, const int * ldb,
			 double * x, const int * ldx, double * rcond, double * ferr, double * berr, double *
			 work, int * iwork, int * info);
  void PREFIX DGETRF_F77(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
  void PREFIX DGEQRF_F77(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info);
  void PREFIX DGETRI_F77(const int* n, double* a, const int* lda, int*ipiv, double * work , const int* lwork, int* info);
  void PREFIX DGETRS_F77(Epetra_fcd, const int* n, const int* nrhs, const double* a, const int* lda, const int* ipiv, double* x ,
			 const int* ldx, int* info);
  void PREFIX DGGEV_F77(Epetra_fcd, Epetra_fcd, const int * n, double * a, const int * lda, double * b, const int * ldb,
			double * alphar, double * alphai, double * beta, double * vl, const int * ldvl,
			double * vr, const int * ldvr, double * work, const int * lwork, int * info);
  void PREFIX DGGLSE_F77(const int * m, const int * n, const int * p, double * a, const int * lda, double * b, const int * ldb,
			 double * c, double * d, double * x, double * work, const int * lwork, int * info);
  void PREFIX DGGSVD_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * m, const int * n, const int * p, int * k, int * l,
			 double * a, const int * lda, double * b, const int * ldb, double * alpha, double * beta,
			 double * u, const int * ldu, double * v, const int * ldv, double * q, const int * ldq, double * work,
    #ifdef HAVE_EPETRA_LAPACK_GSSVD3
			 const int * lwork,
    #endif
			 int * iwork, int * info);
  void PREFIX DHSEQR_F77(Epetra_fcd job, Epetra_fcd, const int * n, const int * ilo, const int * ihi, double * h, const int * ldh,
			 double * wr, double * wi, double * z, const int * ldz, double * work, const int * lwork, int * info);
  double PREFIX DLAMCH_F77(Epetra_fcd);
  void PREFIX DLARFT_F77(Epetra_fcd direct, Epetra_fcd storev, const int * n, const int * k, double * v, const int * ldv, double * tau, double * t, const int * ldt );
  void PREFIX DORGQR_F77(const int * m, const int * n, const int * k, double * a, const int * lda, const double * tau, double * work,
			 const int * lwork, int * info);
  void PREFIX DORGHR_F77(const int * n, const int * ilo, const int * ihi, double * a, const int * lda, const double * tau, double * work,
			 const int * lwork, int * info);
  void PREFIX DORMHR_F77(Epetra_fcd, Epetra_fcd, const int * m, const int * n, const int * ilo, const int * ihi, const double * a,
			 const int * lda, const double * tau, double * c, const int * ldc, double * work, const int * lwork, int * info);
  void PREFIX DPOCON_F77(Epetra_fcd, const int* n, const double* a, const int* lda, const double * anorm, double * rcond,
			 double * work, int * iwork, int* info);
  void PREFIX DPOEQU_F77(const int* n, const double* a, const int* lda, double * s, double * scond, double * amax, int* info);
  void PREFIX DPORFS_F77(Epetra_fcd, const int * n, const int * nrhs, const double * a, const int * lda, const double * af, const int * ldaf,
			 const double * b, const int * ldb, double * x, const int * ldx, double * ferr, double * berr,
			 double * work, int * iwork, int * info);
  void PREFIX DPOSV_F77(Epetra_fcd, const int * n, const int * nrhs, const double* a, const int* lda, double*x , const int* ldx, int* info);
  void PREFIX DPOSVX_F77(Epetra_fcd, Epetra_fcd, const int * n, const int * nrhs, double * a, const int * lda, double * af,
			 const int * ldaf, Epetra_fcd, double * s, double * b, const int * ldb, double * x,
			 const int * ldx, double * rcond, double * ferr, double * berr, double * work,
			 int * iwork, int * info);
  void PREFIX DPOTRF_F77(Epetra_fcd, const int* n, double* a, const int* lda, int* info);
  void PREFIX DPOTRI_F77(Epetra_fcd, const int* n, double* a, const int* lda, int* info);
  void PREFIX DPOTRS_F77(Epetra_fcd, const int * n, const int * nrhs, const double* a, const int* lda, double*x ,
			 const int* ldx, int* info);
  void PREFIX DSPEV_F77( Epetra_fcd, Epetra_fcd,const  int * n, double * ap, double * w, double * z,
			 const int * ldz, double * work, int * info);
  void PREFIX DSPGV_F77(const int * itype, Epetra_fcd, Epetra_fcd, const int * n, double * ap, double * bp,
			double * w, double * z, const int * ldz, double * work, int * info);
  void PREFIX DSTEV_F77(Epetra_fcd jobz, const int * n, double * d, double * e, double * z, const int * ldz,
			double * work, int * info);
  void PREFIX DSYEVD_F77(Epetra_fcd, Epetra_fcd, const int * n, double * a, const int * lda, double * w,
			 double * work, const int * lwork, int * iwork, const int * liwork,int * info);
  void PREFIX DSYEV_F77(Epetra_fcd, Epetra_fcd, const int * n, double * a, const int * lda, double * w,
			double * work, const int * lwork, int * info);
  void PREFIX DSYEVR_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, double * a, const int * lda,
			 const double * vl, const double * vu, const int * il, const int * iu, const
			 double * abstol, int * m,
			 double * w, double * z, const int * ldz,  int * isuppz, double * work,
			 const int * lwork, int * iwork, const int * liwork, int * info);
  void PREFIX DSYEVX_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, double * a, const int * lda,
			 const double * vl, const double * vu, const int * il, const int * iu, const double * abstol, int * m,
			 double * w, double * z, const int * ldz, double * work, const int * lwork, int * iwork,
			 int * ifail, int * info);
  void PREFIX DSYGV_F77(const int * itype, Epetra_fcd, Epetra_fcd, const int * n, double * a, const int * lda,
			double * b, const int * ldb, double * w, double * work, const int * lwork, int * info);
  void PREFIX DSYGVX_F77(const int * itype, Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, double * a,
			 const int * lda, double * b, const int * ldb, const double * vl, const double * vu, const int * il,
			 const int * iu, const double * abstol, int * m, double * w, double * z, const int * ldz,
			 double * work, const int * lwork, int * iwork, int * ifail, int * info);
  void PREFIX DTREVC_F77(Epetra_fcd, Epetra_fcd, int * select, const int * n, const double * t, const int * ldt,
			 double *vl, const int * ldvl, double * vr, const int * ldvr, const int * mm, int * m,
			 double * work, int * info);
  void PREFIX DTREXC_F77(Epetra_fcd, const int * n, double * t, const int * ldt, double * q, const int * ldq,
			 int * ifst, int * ilst, double * work, int * info);
  void PREFIX DTRTRS_F77(Epetra_fcd uplo, Epetra_fcd trans, Epetra_fcd diag, const int *n, const int *nrhs, const double *a,
                         const int *lda, double *b, const int *ldb, int *info);


  void PREFIX SGECON_F77(Epetra_fcd norm, const int* n, const float* a, const int* lda, const float *anorm, float * rcond,
			 float * work, int * iwork, int* info);
  void PREFIX SGEEQU_F77(const int* m, const int* n, const float* a, const int* lda, float * r, float * c, float * rowcnd,
			 float * colcnd, float * amax, int* info);
  void PREFIX SGEEV_F77(Epetra_fcd, Epetra_fcd, const int* n, float* a, const int* lda, float* wr, float* wi,
			float* vl, const int* ldvl,
			float* vr, const int* ldvr, float* work, const int* lwork, int* info);
  void PREFIX SGEEVX_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, float * a, const int * lda,
			 float * wr, float * wi, float * vl, const int * ldvl, float * vr, const int * ldvr,
			 int * ilo, int * ihi, float * scale, float * abnrm, float * rconde, float * rcondv,
			 float * work, const int * lwork, int * iwork, int * info);
  void PREFIX SGEHRD_F77(const int * n, const int * ilo, const int * ihi, float * A, const int * lda, float * tau, float * work,
			 const int * lwork, int * info);
  void PREFIX SGELS_F77(Epetra_fcd ch, const int* m, const int* n, const int* nrhs, float* a, const int* lda, float* b, const int* ldb,
			float* work, const int* lwork, int* info);
  void PREFIX SGELSS_F77(const int * m, const int * n, const int * nrhs, float * a, const int * lda, float * b, const int * ldb,
			 float * s, const float * rcond, int * rank, float * work, const int * lwork, int * info);
  void PREFIX SGEQPF_F77(const int * m, const int * n, float * a, const int * lda, int * jpvt, float * tau, float * work, int * info);
  void PREFIX SGERFS_F77(Epetra_fcd, const int * n, const int * nrhs, const float * a, const int * lda, const float * af, const int * ldaf,
			 const int*ipiv, const float * b, const int * ldb, float * x, const int * ldx, float * ferr, float * berr,
			 float * work, int * iwork, int * info);

  void PREFIX SGESDD_F77(Epetra_fcd, const int * m, const int * n, float * a, const int * lda, float * s, float * u,
			 const int * ldu, float * vt, const int * ldvt, float * work, const int * lwork, int * iwork, int * info);

  void PREFIX SGESVD_F77(Epetra_fcd, Epetra_fcd, const int* m, const int* n, float* a, const int* lda, float* s, float* u,
			 const int* ldu, float* vt, const int* ldvt, float* work, const int* lwork, int* info);
  void PREFIX SGESV_F77(const int * n, const int * nrhs, float* a, const int* lda, int*ipiv, float*x , const int* ldx, int* info);
  void PREFIX SGESVX_F77(Epetra_fcd, Epetra_fcd, const int * n, const int * nrhs, float * a, const int * lda, float * af,
			 const int * ldaf, int*ipiv, Epetra_fcd, float * r, float *c, float * b, const int * ldb,
			 float * x, const int * ldx, float * rcond, float * ferr, float * berr, float *
			 work, int * iwork, int * info);
  void PREFIX SGETRF_F77(const int* m, const int* n, float* a, const int* lda, int* ipiv, int* info);
  void PREFIX SGEQRF_F77(const int* m, const int* n, float* a, const int* lda, float* tau, float* work, const int* lwork, int* info);
  void PREFIX SGETRI_F77(const int* n, float* a, const int* lda, int*ipiv, float * work , const int* lwork, int* info);
  void PREFIX SGETRS_F77(Epetra_fcd, const int* n, const int* nrhs, const float* a, const int* lda, const int* ipiv, float* x ,
			 const int* ldx, int* info);
  void PREFIX SGGEV_F77(Epetra_fcd, Epetra_fcd, const int * n, float * a, const int * lda, float * b, const int * ldb,
			float * alphar, float * alphai, float * beta, float * vl, const int * ldvl,
			float * vr, const int * ldvr, float * work, const int * lwork, int * info);
  void PREFIX SGGLSE_F77(const int * m, const int * n, const int * p, float * a, const int * lda, float * b, const int * ldb,
			 float * c, float * d, float * x, float * work, const int * lwork, int * info);
  void PREFIX SGGSVD_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * m, const int * n, const int * p, int * k, int * l,
			 float * a, const int * lda, float * b, const int * ldb, float * alpha, float * beta,
			 float * u, const int * ldu, float * v, const int * ldv, float * q, const int * ldq, float * work,
    #ifdef HAVE_EPETRA_LAPACK_GSSVD3
			 const int * lwork,
    #endif
			 int * iwork, int * info);
  void PREFIX SHSEQR_F77(Epetra_fcd job, Epetra_fcd, const int * n, const int * ilo, const int * ihi, float * h, const int * ldh,
			 float * wr, float * wi, float * z, const int * ldz, float * work, const int * lwork, int * info);
  float PREFIX SLAMCH_F77(Epetra_fcd);
  void PREFIX SLARFT_F77(Epetra_fcd direct, Epetra_fcd storev, const int * n, const int * k, float * v, const int * ldv, float * tau, float * t, const int * ldt );
  void PREFIX SORGQR_F77(const int * m, const int * n, const int * k, float * a, const int * lda, const float * tau, float * work,
			 const int * lwork, int * info);
  void PREFIX SORGHR_F77(const int * n, const int * ilo, const int * ihi, float * a, const int * lda, const float * tau, float * work,
			 const int * lwork, int * info);
  void PREFIX SORMHR_F77(Epetra_fcd, Epetra_fcd, const int * m, const int * n, const int * ilo, const int * ihi, const float * a,
			 const int * lda, const float * tau, float * c, const int * ldc, float * work, const int * lwork, int * info);
  void PREFIX SPOCON_F77(Epetra_fcd, const int* n, const float* a, const int* lda, const float * anorm, float * rcond,
			 float * work, int * iwork, int* info);
  void PREFIX SPOEQU_F77(const int* n, const float* a, const int* lda, float * s, float * scond, float * amax, int* info);
  void PREFIX SPORFS_F77(Epetra_fcd, const int * n, const int * nrhs, const float * a, const int * lda, const float * af, const int * ldaf,
			 const float * b, const int * ldb, float * x, const int * ldx, float * ferr, float * berr,
			 float * work, int * iwork, int * info);
  void PREFIX SPOSV_F77(Epetra_fcd, const int * n, const int * nrhs, const float* a, const int* lda, float*x , const int* ldx, int* info);
  void PREFIX SPOSVX_F77(Epetra_fcd, Epetra_fcd, const int * n, const int * nrhs, float * a, const int * lda, float * af,
			 const int * ldaf, Epetra_fcd, float * s, float * b, const int * ldb, float * x,
			 const int * ldx, float * rcond, float * ferr, float * berr, float * work,
			 int * iwork, int * info);
  void PREFIX SPOTRF_F77(Epetra_fcd, const int* n, float* a, const int* lda, int* info);
  void PREFIX SPOTRI_F77(Epetra_fcd, const int* n, float* a, const int* lda, int* info);
  void PREFIX SPOTRS_F77(Epetra_fcd, const int * n, const int * nrhs, const float* a, const int* lda, float*x ,
			 const int* ldx, int* info);
  void PREFIX SSPEV_F77( Epetra_fcd, Epetra_fcd,const  int * n, float * ap, float * w, float * z,
			 const int * ldz, float * work, int * info);
  void PREFIX SSPGV_F77(const int * itype, Epetra_fcd, Epetra_fcd, const int * n, float * ap, float * bp,
			float * w, float * z, const int * ldz, float * work, int * info);
  void PREFIX SSTEV_F77(Epetra_fcd jobz, const int * n, float * d, float * e, float * z, const int * ldz,
			float * work, int * info);
  void PREFIX SSYEVD_F77(Epetra_fcd, Epetra_fcd, const int * n, float * a, const int * lda, float * w,
			 float * work, const int * lwork, int * iwork, const int * liwork, int * info);
  void PREFIX SSYEV_F77(Epetra_fcd, Epetra_fcd, const int * n, float * a, const int * lda, float * w,
			float * work, const int * lwork, int * info);
  void PREFIX SSYEVR_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, float * a, const int * lda,
			 const float * vl, const float * vu, const int * il, const int * iu, const
			 float * abstol, int * m,
			 float * w, float * z, const int * ldz,  int * isuppz, float * work,
			 const int * lwork, int * iwork, const int * liwork, int * info);
  void PREFIX SSYEVX_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, float * a, const int * lda,
			 const float * vl, const float * vu, const int * il, const int * iu, const float * abstol, int * m,
			 float * w, float * z, const int * ldz, float * work, const int * lwork, int * iwork,
			 int * ifail, int * info);
  void PREFIX SSYGV_F77(const int * itype, Epetra_fcd, Epetra_fcd, const int * n, float * a, const int * lda,
			float * b, const int * ldb, float * w, float * work, const int * lwork, int * info);
  void PREFIX SSYGVX_F77(const int * itype, Epetra_fcd, Epetra_fcd, Epetra_fcd, const int * n, float * a,
			 const int * lda, float * b, const int * ldb, const float * vl, const float * vu, const int * il,
			 const int * iu, const float * abstol, int * m, float * w, float * z, const int * ldz,
			 float * work, const int * lwork, int * iwork, int * ifail, int * info);
  void PREFIX STREVC_F77(Epetra_fcd, Epetra_fcd, int * select, const int * n, const float * t, const int * ldt,
			 float *vl, const int * ldvl, float * vr, const int * ldvr, const int * mm, int * m,
			 float * work, int * info);
  void PREFIX STREXC_F77(Epetra_fcd, const int * n, float * t, const int * ldt, float * q, const int * ldq,
			 int * ifst, int * ilst, float * work, int * info);
  void PREFIX STRTRS_F77(Epetra_fcd uplo, Epetra_fcd trans, Epetra_fcd diag, const int *n, const int *nrhs, const float *a,
                         const int *lda, float *b, const int *ldb, int *info);


#ifdef __cplusplus
}
#endif

#endif /* EPETRA_LAPACK_WRAPPERS_H */
