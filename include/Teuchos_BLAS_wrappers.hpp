// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_BLAS_WRAPPERS_HPP_
#define _TEUCHOS_BLAS_WRAPPERS_HPP_

#include "Teuchos_ConfigDefs.hpp"
#ifdef _MSC_VER
/* disable warning for C-linkage returning complex class */
#pragma warning ( disable : 4190 )
#endif

/*! \file Teuchos_BLAS_wrappers.hpp
    \brief The Templated BLAS wrappers.
*/


/* A) Define PREFIX and Teuchos_fcd based on platform. */

#if defined(INTEL_CXML)
#  define PREFIX __stdcall
#  define Teuchos_fcd const char *, unsigned int
#elif defined(INTEL_MKL)
#  define PREFIX
#  define Teuchos_fcd const char *
#else /* Not CRAY_T3X or INTEL_CXML or INTEL_MKL */
#  define PREFIX
#  define Teuchos_fcd const char *
#endif

// Handle Complex types (we assume C++11 at a minumum
// Microsoft C++11 apparently does not support float/double _Complex

#if ( defined(_MSC_VER) )
#define Teuchos_Complex_double_type_name std::complex<double>
#define Teuchos_Complex_float_type_name std::complex<float>
#else
#define Teuchos_Complex_double_type_name double _Complex
#define Teuchos_Complex_float_type_name float _Complex
#endif


/* B) Take care of of the link name case */


#define DROTG_F77   F77_BLAS_MANGLE(drotg,DROTG)
#define DROT_F77    F77_BLAS_MANGLE(drot,DROT)
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
#define DSWAP_F77   F77_BLAS_MANGLE(dswap,DSWAP)
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DSYRK_F77   F77_BLAS_MANGLE(dsyrk,DSYRK)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZROTG_F77   F77_BLAS_MANGLE(zrotg,ZROTG)
#define ZROT_F77    F77_BLAS_MANGLE(zrot,ZROT)
#define ZASUM_F77   F77_BLAS_MANGLE(dzasum,DZASUM)
#define ZAXPY_F77   F77_BLAS_MANGLE(zaxpy,ZAXPY)
#define ZCOPY_F77   F77_BLAS_MANGLE(zcopy,ZCOPY)
#define ZDOT_F77    F77_BLAS_MANGLE(zdotc,ZDOTC)
#define ZNRM2_F77   F77_BLAS_MANGLE(dznrm2,DZNRM2)
#define ZSCAL_F77   F77_BLAS_MANGLE(zscal,ZSCAL)
#define IZAMAX_F77  F77_BLAS_MANGLE(izamax,IZAMAX)
#define ZGEMV_F77   F77_BLAS_MANGLE(zgemv,ZGEMV)
#define ZGER_F77    F77_BLAS_MANGLE(zgeru,ZGERU)
#define ZTRMV_F77   F77_BLAS_MANGLE(ztrmv,ZTRMV)
#define ZGEMM_F77   F77_BLAS_MANGLE(zgemm,ZGEMM)
#define ZSWAP_F77   F77_BLAS_MANGLE(zswap,ZSWAP)
#define ZSYMM_F77   F77_BLAS_MANGLE(zsymm,ZSYMM)
#define ZSYRK_F77   F77_BLAS_MANGLE(zsyrk,ZSYRK)
#define ZHERK_F77   F77_BLAS_MANGLE(zherk,ZHERK)
#define ZTRMM_F77   F77_BLAS_MANGLE(ztrmm,ZTRMM)
#define ZTRSM_F77   F77_BLAS_MANGLE(ztrsm,ZTRSM)

#endif /* HAVE_TEUCHOS_COMPLEX */

#define SROTG_F77   F77_BLAS_MANGLE(srotg,SROTG)
#define SROT_F77    F77_BLAS_MANGLE(srot,SROT)
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
#define SSWAP_F77   F77_BLAS_MANGLE(sswap,SSWAP)
#define SSYMM_F77   F77_BLAS_MANGLE(ssymm,SSYMM)
#define SSYRK_F77   F77_BLAS_MANGLE(ssyrk,SSYRK)
#define STRMM_F77   F77_BLAS_MANGLE(strmm,STRMM)
#define STRSM_F77   F77_BLAS_MANGLE(strsm,STRSM)

#ifdef HAVE_TEUCHOS_COMPLEX

#define CROTG_F77   F77_BLAS_MANGLE(crotg,CROTG)
#define CROT_F77    F77_BLAS_MANGLE(crot,CROT)
#define SCASUM_F77  F77_BLAS_MANGLE(scasum,SCASUM)
#define CAXPY_F77   F77_BLAS_MANGLE(caxpy,CAXPY)
#define CCOPY_F77   F77_BLAS_MANGLE(ccopy,CCOPY)
#define CDOT_F77    F77_BLAS_MANGLE(cdotc,CDOTC)
#define SCNRM2_F77   F77_BLAS_MANGLE(scnrm2,SCNRM2)
#define CSCAL_F77   F77_BLAS_MANGLE(cscal,CSCAL)
#define ICAMAX_F77  F77_BLAS_MANGLE(icamax,ICAMAX)
#define CGEMV_F77   F77_BLAS_MANGLE(cgemv,CGEMV)
#define CGER_F77    F77_BLAS_MANGLE(cgeru,CGERU)
#define CTRMV_F77   F77_BLAS_MANGLE(ctrmv,CTRMV)
#define CGEMM_F77   F77_BLAS_MANGLE(cgemm,CGEMM)
#define CSWAP_F77   F77_BLAS_MANGLE(cswap,CSWAP)
#define CSYMM_F77   F77_BLAS_MANGLE(csymm,CSYMM)
#define CSYRK_F77   F77_BLAS_MANGLE(csyrk,CSYRK)
#define CHERK_F77   F77_BLAS_MANGLE(cherk,CHERK)
#define CTRMM_F77   F77_BLAS_MANGLE(ctrmm,CTRMM)
#define CTRSM_F77   F77_BLAS_MANGLE(ctrsm,CTRSM)
#define TEUCHOS_BLAS_CONVERT_COMPLEX_FORTRAN_TO_CXX(TYPE, Z) \
   reinterpret_cast<std::complex<TYPE>&>(Z);
// NOTE: The above is guaranteed to be okay given the C99 and C++11 standards

#endif /* HAVE_TEUCHOS_COMPLEX */


/* C) Define the function prototypes for all platforms! */

#ifdef __cplusplus
extern "C" {
#endif


/* Double precision BLAS 1 */
void PREFIX DROTG_F77(double* da, double* db, double* c, double* s);
void PREFIX DROT_F77(const int* n, double* dx, const int* incx, double* dy, const int* incy, double* c, double* s);
double PREFIX DASUM_F77(const int* n, const double x[], const int* incx);
void PREFIX DAXPY_F77(const int* n, const double* alpha, const double x[], const int* incx, double y[], const int* incy);
void PREFIX DCOPY_F77(const int* n, const double *x, const int* incx, double *y, const int* incy);
double PREFIX DDOT_F77(const int* n, const double x[], const int* incx, const double y[], const int* incy);
double PREFIX DNRM2_F77(const int* n, const double x[], const int* incx);
void PREFIX DSCAL_F77(const int* n, const double* alpha, double *x, const int* incx);
void PREFIX DSWAP_F77(const int* const n, double* const x, const int* const incx,
                      double* const y, const int* const incy);
int PREFIX IDAMAX_F77(const int* n, const double *x, const int* incx);

/* Double std::complex precision BLAS 1 */
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus)

#  if defined(HAVE_COMPLEX_BLAS_PROBLEM)
#    if defined(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
void PREFIX ZDOT_F77(std::complex<double> *ret, const int* n, const std::complex<double> x[], const int* incx, const std::complex<double> y[], const int* incy);
#    elif defined(HAVE_VECLIB_COMPLEX_BLAS)
// no declarations; they're in cblas.h
#      include <vecLib/cblas.h>
#    else
  // mfh 01 Feb 2013: If the code reaches this point, it means that
  // some complex BLAS routines are broken, but there is no easy
  // workaround.  We deal with this in Teuchos_BLAS.cpp by
  // reimplementing the offending routines.
#    endif // HAVE_COMPLEX_BLAS_PROBLEM
#  else // no problem
Teuchos_Complex_double_type_name PREFIX ZDOT_F77(const int* n, const std::complex<double> x[], const int* incx, const std::complex<double> y[], const int* incy);
#  endif // defined(HAVE_COMPLEX_BLAS_PROBLEM)

double PREFIX ZNRM2_F77(const int* n, const std::complex<double> x[], const int* incx);
double PREFIX ZASUM_F77(const int* n, const std::complex<double> x[], const int* incx);
void PREFIX ZROTG_F77(std::complex<double>* da, std::complex<double>* db, double* c, std::complex<double>* s);
void PREFIX ZROT_F77(const int* n, std::complex<double>* dx, const int* incx, std::complex<double>* dy, const int* incy, double* c, std::complex<double>* s);
void PREFIX ZAXPY_F77(const int* n, const std::complex<double>* alpha, const std::complex<double> x[], const int* incx, std::complex<double> y[], const int* incy);
void PREFIX ZCOPY_F77(const int* n, const std::complex<double> *x, const int* incx, std::complex<double> *y, const int* incy);
void PREFIX ZSCAL_F77(const int* n, const std::complex<double>* alpha, std::complex<double> *x, const int* incx);
void PREFIX ZSWAP_F77(const int* const n, std::complex<double>* const x, const int* const incx,
                      std::complex<double>* const y, const int* const incy);
int PREFIX IZAMAX_F77(const int* n, const std::complex<double> *x, const int* incx);

#endif /* defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus) */

/* Single precision BLAS 1 */
#ifdef HAVE_TEUCHOS_BLASFLOAT
#  ifdef HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX
#    include <vecLib/cblas.h>
#  elif defined(HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN)
double PREFIX SASUM_F77(const int* n, const float x[], const int* incx);
double PREFIX SDOT_F77(const int* n, const float x[], const int* incx, const float y[], const int* incy);
double PREFIX SNRM2_F77(const int* n, const float x[], const int* incx);
#  else
float PREFIX SASUM_F77(const int* n, const float x[], const int* incx);
float PREFIX SDOT_F77(const int* n, const float x[], const int* incx, const float y[], const int* incy);
float PREFIX SNRM2_F77(const int* n, const float x[], const int* incx);
#  endif // which blasfloat
#endif // ifdef blasfloat
void PREFIX SROTG_F77(float* da, float* db, float* c, float* s);
void PREFIX SROT_F77(const int* n, float* dx, const int* incx, float* dy, const int* incy, float* c, float* s);
void PREFIX SAXPY_F77(const int* n, const float* alpha, const float x[], const int* incx, float y[], const int* incy);
void PREFIX SCOPY_F77(const int* n, const float *x, const int* incx, float *y, const int* incy);
void PREFIX SSCAL_F77(const int* n, const float* alpha, float *x, const int* incx);
void PREFIX SSWAP_F77(const int* const n, float* const x, const int* const incx,
                      float* const y, const int* const incy);
int PREFIX ISAMAX_F77(const int* n, const float *x, const int* incx);

/* Single std::complex precision BLAS 1 */
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus)
#  if defined(HAVE_TEUCHOS_BLASFLOAT)
#    if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
// no declarations; they're in cblas.h
#      include <vecLib/cblas.h>
#    elif defined(HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN)
double PREFIX SCASUM_F77(const int* n, const std::complex<float> x[], const int* incx);
double PREFIX SCNRM2_F77(const int* n, const std::complex<float> x[], const int* incx);
#    else
float PREFIX SCASUM_F77(const int* n, const std::complex<float> x[], const int* incx);
float PREFIX SCNRM2_F77(const int* n, const std::complex<float> x[], const int* incx);
#  endif // Whether or not we have the veclib bugfix
#endif // defined(HAVE_TEUCHOS_BLASFLOAT)

#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
// no declarations; they're in cblas.h
#include <vecLib/cblas.h>
#elif defined(HAVE_COMPLEX_BLAS_PROBLEM) && defined(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
void PREFIX CDOT_F77(std::complex<float> *ret, const int* n, const std::complex<float> x[], const int* incx, const std::complex<float> y[], const int* incy);
#elif defined(HAVE_TEUCHOS_BLASFLOAT)
Teuchos_Complex_float_type_name PREFIX CDOT_F77(const int* n, const std::complex<float> x[], const int* incx, const std::complex<float> y[], const int* incy);
#else
// the code is literally in Teuchos_BLAS.cpp
#endif

void PREFIX CROTG_F77(std::complex<float>* da, std::complex<float>* db, float* c, std::complex<float>* s);
void PREFIX CROT_F77(const int* n, std::complex<float>* dx, const int* incx, std::complex<float>* dy, const int* incy, float* c, std::complex<float>* s);
void PREFIX CAXPY_F77(const int* n, const std::complex<float>* alpha, const std::complex<float> x[], const int* incx, std::complex<float> y[], const int* incy);
void PREFIX CCOPY_F77(const int* n, const std::complex<float> *x, const int* incx, std::complex<float> *y, const int* incy);
void PREFIX CSCAL_F77(const int* n, const std::complex<float>* alpha, std::complex<float> *x, const int* incx);
void PREFIX CSWAP_F77(const int* const n, std::complex<float>* const x, const int* const incx,
                      std::complex<float>* const y, const int* const incy);
int PREFIX ICAMAX_F77(const int* n, const std::complex<float> *x, const int* incx);

#endif /* defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus) */

/* Double precision BLAS 2 */
void PREFIX DGEMV_F77(Teuchos_fcd, const int* m, const int* n, const double* alpha, const double A[], const int* lda,
                 const double x[], const int* incx, const double* beta, double y[], const int* incy);
void PREFIX DTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n,
                const double *a, const int *lda, double *x, const int *incx);
void PREFIX DGER_F77(const int *m, const int *n, const double *alpha, const double *x, const int *incx, const double *y,
               const int *incy, double *a, const int *lda);

/* Double precision BLAS 2 */
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus)

void PREFIX ZGEMV_F77(Teuchos_fcd, const int* m, const int* n, const std::complex<double>* alpha, const std::complex<double> A[], const int* lda,
                 const std::complex<double> x[], const int* incx, const std::complex<double>* beta, std::complex<double> y[], const int* incy);
void PREFIX ZTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n,
                const std::complex<double> *a, const int *lda, std::complex<double> *x, const int *incx);
void PREFIX ZGER_F77(const int *m, const int *n, const std::complex<double> *alpha, const std::complex<double> *x, const int *incx, const std::complex<double> *y,
               const int *incy, std::complex<double> *a, const int *lda);

#endif /* defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus) */

/* Single precision BLAS 2 */
void PREFIX SGEMV_F77(Teuchos_fcd, const int* m, const int* n, const float* alpha, const float A[], const int* lda,
                 const float x[], const int* incx, const float* beta, float y[], const int* incy);
void PREFIX STRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n,
                const float *a, const int *lda, float *x, const int *incx);
void PREFIX SGER_F77(const int *m, const int *n, const float *alpha, const float *x, const int *incx, const float *y,
               const int *incy, float *a, const int *lda);

/* Single std::complex precision BLAS 2 */
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus)

void PREFIX CGEMV_F77(Teuchos_fcd, const int* m, const int* n, const std::complex<float>* alpha, const std::complex<float> A[], const int* lda,
                 const std::complex<float> x[], const int* incx, const std::complex<float>* beta, std::complex<float> y[], const int* incy);
void PREFIX CTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n,
                const std::complex<float> *a, const int *lda, std::complex<float> *x, const int *incx);
void PREFIX CGER_F77(const int *m, const int *n, const std::complex<float> *alpha, const std::complex<float> *x, const int *incx, const std::complex<float> *y,
               const int *incy, std::complex<float> *a, const int *lda);

#endif /* defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus) */

/* Double precision BLAS 3 */
void PREFIX DGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int *
                n, const int *k, const double *alpha, const double *a, const int *lda,
                const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
void PREFIX DSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const double *alpha, const double *a, const int *lda,
                const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
void PREFIX DSYRK_F77(Teuchos_fcd, Teuchos_fcd, const int *n, const int * k,
                const double *alpha, const double *a, const int *lda,
                const double *beta, double *c, const int *ldc);
void PREFIX DTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const double *alpha, const double *a, const int * lda, double *b, const int *ldb);
void PREFIX DTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const double *alpha, const double *a, const int *
                lda, double *b, const int *ldb);

/* Double std::complex precision BLAS 3 */
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus)

void PREFIX ZGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int *
                n, const int *k, const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
                const std::complex<double> *b, const int *ldb, const std::complex<double> *beta, std::complex<double> *c, const int *ldc);
void PREFIX ZSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
                const std::complex<double> *b, const int *ldb, const std::complex<double> *beta, std::complex<double> *c, const int *ldc);
void PREFIX ZSYRK_F77(Teuchos_fcd, Teuchos_fcd, const int *n, const int * k,
                const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
                const std::complex<double> *beta, std::complex<double> *c, const int *ldc);
void PREFIX ZHERK_F77(Teuchos_fcd, Teuchos_fcd, const int *n, const int * k,
                const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
                const std::complex<double> *beta, std::complex<double> *c, const int *ldc);
void PREFIX ZTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const std::complex<double> *alpha, const std::complex<double> *a, const int * lda, std::complex<double> *b, const int *ldb);
void PREFIX ZTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const std::complex<double> *alpha, const std::complex<double> *a, const int *
                lda, std::complex<double> *b, const int *ldb);

#endif /* defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus) */

/* Single precision BLAS 3 */
void PREFIX SGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int *
                n, const int *k, const float *alpha, const float *a, const int *lda,
                const float *b, const int *ldb, const float *beta, float *c, const int *ldc);
void PREFIX SSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const float *alpha, const float *a, const int *lda,
                const float *b, const int *ldb, const float *beta, float *c, const int *ldc);
void PREFIX SSYRK_F77(Teuchos_fcd, Teuchos_fcd, const int *n, const int * k,
                const float *alpha, const float *a, const int *lda,
                const float *beta, float *c, const int *ldc);
void PREFIX STRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const float *alpha, const float *a, const int * lda, float *b, const int *ldb);
void PREFIX STRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const float *alpha, const float *a, const int *
                lda, float *b, const int *ldb);

/* Single std::complex precision BLAS 3 */

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus)

void PREFIX CGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int *
                n, const int *k, const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
                const std::complex<float> *b, const int *ldb, const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
void PREFIX CSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n,
                const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
                const std::complex<float> *b, const int *ldb, const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
void PREFIX CTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const std::complex<float> *alpha, const std::complex<float> *a, const int * lda, std::complex<float> *b, const int *ldb);
void PREFIX CSYRK_F77(Teuchos_fcd, Teuchos_fcd, const int *n, const int * k,
                const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
                const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
void PREFIX CHERK_F77(Teuchos_fcd, Teuchos_fcd, const int *n, const int * k,
                const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
                const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
void PREFIX CTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                const int *m, const int *n, const std::complex<float> *alpha, const std::complex<float> *a, const int *
                lda, std::complex<float> *b, const int *ldb);

#endif /* defined(HAVE_TEUCHOS_COMPLEX) && defined(__cplusplus) */

#ifdef __cplusplus
}
#endif

/* Don't leave a global macros called PREFIX or Teuchos_fcd laying around */

#ifdef PREFIX
#undef PREFIX
#endif

#ifdef Teuchos_fcd
#undef Teuchos_fcd
#endif

#endif /* end of TEUCHOS_BLAS_WRAPPERS_HPP_ */
