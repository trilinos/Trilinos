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

#ifndef EPETRA_BLAS_WRAPPERS_H
#define EPETRA_BLAS_WRAPPERS_H

#include "Epetra_ConfigDefs.h"
/* #include <stdio.h> */
/* #include <string.h> */


/* Define fcd (Fortran Epetra_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)


#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define Epetra_fcd fcd

#define DASUM_F77   SASUM
#define DAXPY_F77   SAXPY
#define DCOPY_F77   SCOPY
#define DDOT_F77    SDOT
#define DNRM2_F77   SNRM2
#define DSCAL_F77   SSCAL
#define IDAMAX_F77  ISAMAX
#define DGEMV_F77   SGEMV
#define DGER_F77    SGER
#define DTRMV_F77   STRMV
#define DGEMM_F77   SGEMM
#define DSYMM_F77   SSYMM
#define DTRMM_F77   STRMM
#define DTRSM_F77   STRSM
#define DSYRK_F77   SSYRK
#define EPETRA_DCRSMV_F77   EPETRA_DCRSMV
#define EPETRA_DCRSMM_F77   EPETRA_DCRSMM
#define EPETRA_DCRSSV_F77   EPETRA_DCRSSV
#define EPETRA_DCRSSM_F77   EPETRA_DCRSSM

#elif defined(INTEL_CXML)

#define PREFIX __stdcall
#define Epetra_fcd const char *, const unsigned int

#elif defined(INTEL_MKL)

#define PREFIX
#define Epetra_fcd const char *

#endif 

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_BLAS_MANGLE is defined undefine it because we want to redefine */

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef F77_FUNC_
#undef F77_FUNC_
#endif

#define F77_BLAS_MANGLE(lcase,UCASE) UCASE
#define F77_FUNC(lcase,UCASE) UCASE
#define F77_FUNC_(lcase,UCASE) UCASE

#else /* Define Epetra_fcd for all other machines */

#define PREFIX
#define Epetra_fcd const char * 

/* Use autoconf's definition of F77_BLAS_MANGLE 
   unless using old make system */

#ifdef TRILINOS_NO_CONFIG_H

#ifdef F77_BLAS_MANGLE
#undef F77_BLAS_MANGLE
#endif
#ifdef F77_FUNC
#undef F77_FUNC
#endif
#ifdef F77_FUNC_
#undef F77_FUNC_
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_BLAS_MANGLE(lcase,UCASE) lcase
#define F77_FUNC(lcase,UCASE) lcase
#define F77_FUNC_(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_BLAS_MANGLE(lcase,UCASE) lcase ## _
#define F77_FUNC(lcase,UCASE) lcase ## _
#define F77_FUNC_(lcase,UCASE) lcase ## __
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */

#endif /* TRILINOS_NO_CONFIG_H */

#endif /* defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL) */

#ifndef CRAY_T3X /* Double declarations already done for the Cray */

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
#define DSYRK_F77   F77_BLAS_MANGLE(dsyrk,DSYRK)

#ifndef FORTRAN_DISABLED

#if defined(__GNUC__) || defined(_WIN32) /* standard Epetra implementation */ 

#define EPETRA_DCRSMV_F77   F77_FUNC_(epetra_dcrsmv,EPETRA_DCRSMV)
#define EPETRA_DCRSMM_F77   F77_FUNC_(epetra_dcrsmm,EPETRA_DCRSMM)
#define EPETRA_DCRSSV_F77   F77_FUNC_(epetra_dcrssv,EPETRA_DCRSSV)
#define EPETRA_DCRSSM_F77   F77_FUNC_(epetra_dcrssm,EPETRA_DCRSSM)

#else /* MSE: 3/17/05 - patch for Solaris/OSF/IRIX */ 

#define EPETRA_DCRSMV_F77   F77_FUNC(epetra_dcrsmv,EPETRA_DCRSMV) 
#define EPETRA_DCRSMM_F77   F77_FUNC(epetra_dcrsmm,EPETRA_DCRSMM) 
#define EPETRA_DCRSSV_F77   F77_FUNC(epetra_dcrssv,EPETRA_DCRSSV) 
#define EPETRA_DCRSSM_F77   F77_FUNC(epetra_dcrssm,EPETRA_DCRSSM) 
#endif /* __GNUC__ */ 

#endif /* FORTRAN_DISABLED */

/* End of defines for double precision when not on a T3X */

#endif

/* The following defines are good for all platforms */


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
#define SSYRK_F77   F77_BLAS_MANGLE(ssyrk,SSYRK)
    
/* Explicitly define each F77 name for all BLAS kernels */

#ifdef __cplusplus
extern "C" {
#endif

/* Double precision BLAS 1 */
double PREFIX DASUM_F77(const int* n, const double x[], const int* incx);
void PREFIX DAXPY_F77(const int* n, const double* alpha, const double x[], const int* incx, double y[], const int* incy);
void PREFIX DCOPY_F77(const int* n, const double *x, const int* incx, double *y, const int* incy);
double PREFIX DDOT_F77(const int* n, const double x[], const int* incx, const double y[], const int* incy);
double PREFIX DNRM2_F77(const int* n, const double x[], const int* incx);
void PREFIX DSCAL_F77(const int* n, const double* alpha, double *x, const int* incx);
int PREFIX IDAMAX_F77(const int* n, const double *x, const int* incx);

/* Single precision BLAS 1 */
float PREFIX SASUM_F77(const int* n, const float x[], const int* incx);
void PREFIX SAXPY_F77(const int* n, const float* alpha, const float x[], const int* incx, float y[], const int* incy);
void PREFIX SCOPY_F77(const int* n, const float *x, const int* incx, float *y, const int* incy);
float PREFIX SDOT_F77(const int* n, const float x[], const int* incx, const float y[], const int* incy);
float PREFIX SNRM2_F77(const int* n, const float x[], const int* incx);
void PREFIX SSCAL_F77(const int* n, const float* alpha, float *x, const int* incx);
int PREFIX ISAMAX_F77(const int* n, const float *x, const int* incx);

/* Double precision BLAS 2 */
void PREFIX DGEMV_F77(Epetra_fcd, const int* m, const int* n, const double* alpha, const double A[], const int* lda,
		       const double x[], const int* incx, const double* beta, double y[], const int* incy);
void PREFIX DTRMV_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *n, 
		      const double *a, const int *lda, double *x, const int *incx);
void PREFIX DGER_F77(const int *m, const int *n, const double *alpha, const double *x, const int *incx, const double *y, 
		     const int *incy, double *a, const int *lda);


/* Single precision BLAS 2 */
void PREFIX SGEMV_F77(Epetra_fcd, const int* m, const int* n, const float* alpha, const float A[], const int* lda,
		       const float x[], const int* incx, const float* beta, float y[], const int* incy);
void PREFIX STRMV_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *n, 
		      const float *a, const int *lda, float *x, const int *incx);
void PREFIX SGER_F77(const int *m, const int *n, const float *alpha, const float *x, const int *incx, const float *y, 
		     const int *incy, float *a, const int *lda);

/* Double precision BLAS 3 */
void PREFIX DGEMM_F77(Epetra_fcd, Epetra_fcd, const int *m, const int *
		      n, const int *k, const double *alpha, const double *a, const int *lda, 
		      const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
void PREFIX DSYMM_F77(Epetra_fcd, Epetra_fcd, const int *m, const int * n,
		      const double *alpha, const double *a, const int *lda, 
		      const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
void PREFIX DTRMM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      const int *m, const int *n, const double *alpha, const double *a, const int * lda, double *b, const int *ldb);
void PREFIX DTRSM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      const int *m, const int *n, const double *alpha, const double *a, const int *
		      lda, double *b, const int *ldb);
void PREFIX EPETRA_DCRSMV_F77(const int *, const int *, const int *, const double *, const int *, 
			      const int *, double *, double *);
void PREFIX EPETRA_DCRSMM_F77(const int *, const int *, const int *, const double *, const int *, 
			      const int *, double *, int *, double *, int *, int *);
void PREFIX EPETRA_DCRSSV_F77(const int *, const int *, const int *, const int *, const int *, 
			      const int *, const double *, const int *, const int *, double *, 
			      double *, const int *);
void PREFIX EPETRA_DCRSSM_F77(const int *, const int *, const int *, const int *, const int *, 
			      const int *, const double *, const int *, const int *, double *, 
			      const int *, double *, const int *, const int *, const int *);
void PREFIX DSYRK_F77(Epetra_fcd uplo, Epetra_fcd trans, const int *n, const int *k,
                  const double *alpha, const double *a, const int *lda, const double *beta,
                  double *c, const int *ldc);

/* Single precision BLAS 3 */
void PREFIX SGEMM_F77(Epetra_fcd, Epetra_fcd, const int *m, const int *
		      n, const int *k, const float *alpha, const float *a, const int *lda, 
		      const float *b, const int *ldb, const float *beta, float *c, const int *ldc);
void PREFIX SSYMM_F77(Epetra_fcd, Epetra_fcd, const int *m, const int * n,
		      const float *alpha, const float *a, const int *lda, 
		      const float *b, const int *ldb, const float *beta, float *c, const int *ldc);
void PREFIX STRMM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      const int *m, const int *n, const float *alpha, const float *a, const int * lda, float *b, const int *ldb);
void PREFIX STRSM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      const int *m, const int *n, const float *alpha, const float *a, const int *
		      lda, float *b, const int *ldb);

void PREFIX XERBLA_F77(Epetra_fcd, int *info);

void PREFIX SSYRK_F77(Epetra_fcd uplo, Epetra_fcd trans, const int *n, const int *k,
                      const float *alpha, const float *a, const int *lda, const float *beta,
                      float *c, const int *ldc);

#ifdef __cplusplus
}
#endif

#endif /* EPETRA_BLAS_WRAPPERS_H */
