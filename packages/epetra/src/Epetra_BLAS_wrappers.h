
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

#ifndef NO_AUTOTOOLS
#include "Epetra_config.h"
#endif /*This file doesn't exist in old make and configure system*/
#include <stdio.h>
#include <string.h>


/* Define fcd (Fortran Epetra_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)


#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define Epetra_fcd fcd

#define DASUM_F77   F77_FUNC(sasum,SASUM)
#define DAXPY_F77   F77_FUNC(saxpy,SAXPY)
#define DCOPY_F77   F77_FUNC(scopy,SCOPY)
#define DDOT_F77    F77_FUNC(sdot,SDOT)
#define DNRM2_F77   F77_FUNC(snrm2,SNRM2)
#define DSCAL_F77   F77_FUNC(sscal,SSCAL)
#define IDAMAX_F77  F77_FUNC(isamax,ISAMAX)
#define DGEMV_F77   F77_FUNC(sgemv,SGEMV)
#define DGER_F77    F77_FUNC(sger,SGER)
#define DTRMV_F77   F77_FUNC(strmv,STRMV)
#define DGEMM_F77   F77_FUNC(sgemm,SGEMM)
#define DSYMM_F77   F77_FUNC(ssymm,SSYMM)
#define DTRMM_F77   F77_FUNC(strmm,STRMM)
#define DTRSM_F77   F77_FUNC(strsm,STRSM)

#elif defined(INTEL_CXML)

#define PREFIX __stdcall
#define Epetra_fcd char *, unsigned int

#define DASUM_F77   F77_FUNC(dasum,DASUM)
#define DAXPY_F77   F77_FUNC(daxpy,DAXPY)
#define DCOPY_F77   F77_FUNC(dcopy,DCOPY)
#define DDOT_F77    F77_FUNC(ddot,DDOT)
#define DNRM2_F77   F77_FUNC(dnrm2,DNRM2)
#define DSCAL_F77   F77_FUNC(dscal,DSCAL)
#define IDAMAX_F77  F77_FUNC(idamax,IDAMAX)
#define DGEMV_F77   F77_FUNC(dgemv,DGEMV)
#define DGER_F77    F77_FUNC(dger,DGER)
#define DTRMV_F77   F77_FUNC(dtrmv,DTRMV)
#define DGEMM_F77   F77_FUNC(dgemm,DGEMM)
#define DSYMM_F77   F77_FUNC(dsymm,DSYMM)
#define DTRMM_F77   F77_FUNC(dtrmm,DTRMM)
#define DTRSM_F77   F77_FUNC(dtrsm,DTRSM)


#elif defined(INTEL_MKL)

#define PREFIX
#define Epetra_fcd char *

#define DASUM_F77   F77_FUNC(dasum,DASUM)
#define DAXPY_F77   F77_FUNC(daxpy,DAXPY)
#define DCOPY_F77   F77_FUNC(dcopy,DCOPY)
#define DDOT_F77    F77_FUNC(ddot,DDOT)
#define DNRM2_F77   F77_FUNC(dnrm2,DNRM2)
#define DSCAL_F77   F77_FUNC(dscal,DSCAL)
#define IDAMAX_F77  F77_FUNC(idamax,IDAMAX)
#define DGEMV_F77   F77_FUNC(dgemv,DGEMV)
#define DGER_F77    F77_FUNC(dger,DGER)
#define DTRMV_F77   F77_FUNC(dtrmv,DTRMV)
#define DGEMM_F77   F77_FUNC(dgemm,DGEMM)
#define DSYMM_F77   F77_FUNC(dsymm,DSYMM)
#define DTRMM_F77   F77_FUNC(dtrmm,DTRMM)
#define DTRSM_F77   F77_FUNC(dtrsm,DTRSM)


#endif 

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_FUNC is defined undefine it because we want to redefine */

#ifdef F77_FUNC
#undef F77_FUNC
#endif


#define F77_FUNC(lcase,UCASE) PREFIX UCASE

#else /* Define Epetra_fcd for all other machines */

#define PREFIX
#define Epetra_fcd char * 

/* Use autoconf's definition of F77_FUNC 
   unless using old make system */

#ifdef NO_AUTOTOOLS

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_FUNC(lcase,UCASE) PREFIX lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_FUNC(lcase,UCASE) PREFIX lcase ## _
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */
#endif /* NO_AUTOTOOLS */

#define DASUM_F77   F77_FUNC(dasum,DASUM)
#define DAXPY_F77   F77_FUNC(daxpy,DAXPY)
#define DCOPY_F77   F77_FUNC(dcopy,DCOPY)
#define DDOT_F77    F77_FUNC(ddot,DDOT)
#define DNRM2_F77   F77_FUNC(dnrm2,DNRM2)
#define DSCAL_F77   F77_FUNC(dscal,DSCAL)
#define IDAMAX_F77  F77_FUNC(idamax,IDAMAX)
#define DGEMV_F77   F77_FUNC(dgemv,DGEMV)
#define DGER_F77    F77_FUNC(dger,DGER)
#define DTRMV_F77   F77_FUNC(dtrmv,DTRMV)
#define DGEMM_F77   F77_FUNC(dgemm,DGEMM)
#define DSYMM_F77   F77_FUNC(dsymm,DSYMM)
#define DTRMM_F77   F77_FUNC(dtrmm,DTRMM)
#define DTRSM_F77   F77_FUNC(dtrsm,DTRSM)


#endif


#define SSCAL_F77   F77_FUNC(sscal,SSCAL)
#define SCOPY_F77   F77_FUNC(scopy,SCOPY)
#define SAXPY_F77   F77_FUNC(saxpy,SAXPY)
#define SDOT_F77    F77_FUNC(sdot,SDOT)
#define SNRM2_F77   F77_FUNC(snrm2,SNRM2)
#define SASUM_F77   F77_FUNC(sasum,SASUM)
#define ISAMAX_F77  F77_FUNC(isamax,ISAMAX)

#define SGEMV_F77   F77_FUNC(sgemv,SGEMV)
#define SGER_F77    F77_FUNC(sger,SGER)
#define STRMV_F77   F77_FUNC(strmv,STRMV)
#define SGEMM_F77   F77_FUNC(sgemm,SGEMM)
#define SSYMM_F77   F77_FUNC(ssymm,SSYMM)
#define STRMM_F77   F77_FUNC(strmm,STRMM)
#define STRSM_F77   F77_FUNC(strsm,STRSM)
    
/* Explicitly define each F77 name for all BLAS kernels */

#ifdef __cplusplus
extern "C" {
#endif

/* Double precision BLAS 1 */
double DASUM_F77(int* n, double x[], int* incx);
void DAXPY_F77(int* n, double* alpha, double x[], int* incx, double y[], int* incy);
void DCOPY_F77(int* n, double *x, int* incx, double *y, int* incy);
double DDOT_F77(int* n, double x[], int* incx, double y[], int* incy);
double DNRM2_F77(int* n, double x[], int* incx);
void DSCAL_F77(int* n, double* alpha, double *x, int* incx);
int IDAMAX_F77(int* n, double *x, int* incx);

/* Single precision BLAS 1 */
float SASUM_F77(int* n, float x[], int* incx);
void SAXPY_F77(int* n, float* alpha, float x[], int* incx, float y[], int* incy);
void SCOPY_F77(int* n, float *x, int* incx, float *y, int* incy);
float SDOT_F77(int* n, float x[], int* incx, float y[], int* incy);
float SNRM2_F77(int* n, float x[], int* incx);
void SSCAL_F77(int* n, float* alpha, float *x, int* incx);
int ISAMAX_F77(int* n, float *x, int* incx);

/* Double precision BLAS 2 */
void DGEMV_F77(Epetra_fcd, int* m, int* n, double* alpha, double A[], int* lda,
		       double x[], int* incx, double* beta, double y[], int* incy);
void DTRMV_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, int *n, 
		      double *a, int *lda, double *x, int *incx);
void DGER_F77(int *m, int *n, double *alpha, double *x, int *incx, double *y, 
		     int *incy, double *a, int *lda);


/* Single precision BLAS 2 */
void SGEMV_F77(Epetra_fcd, int* m, int* n, float* alpha, float A[], int* lda,
		       float x[], int* incx, float* beta, float y[], int* incy);
void STRMV_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, int *n, 
		      float *a, int *lda, float *x, int *incx);
void SGER_F77(int *m, int *n, float *alpha, float *x, int *incx, float *y, 
		     int *incy, float *a, int *lda);

/* Double precision BLAS 3 */
void DGEMM_F77(Epetra_fcd, Epetra_fcd, int *m, int *
		      n, int *k, double *alpha, double *a, int *lda, 
		      double *b, int *ldb, double *beta, double *c, int *ldc);
void DSYMM_F77(Epetra_fcd, Epetra_fcd, int *m, int * n,
		      double *alpha, double *a, int *lda, 
		      double *b, int *ldb, double *beta, double *c, int *ldc);
void DTRMM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, double *alpha, double *a, int * lda, double *b, int *ldb);
void DTRSM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, double *alpha, double *a, int *
		      lda, double *b, int *ldb);

/* Single precision BLAS 3 */
void SGEMM_F77(Epetra_fcd, Epetra_fcd, int *m, int *
		      n, int *k, float *alpha, float *a, int *lda, 
		      float *b, int *ldb, float *beta, float *c, int *ldc);
void SSYMM_F77(Epetra_fcd, Epetra_fcd, int *m, int * n,
		      float *alpha, float *a, int *lda, 
		      float *b, int *ldb, float *beta, float *c, int *ldc);
void STRMM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, float *alpha, float *a, int * lda, float *b, int *ldb);
void STRSM_F77(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, float *alpha, float *a, int *
		      lda, float *b, int *ldb);

void XERBLA_F77(Epetra_fcd, int *info);

#ifdef __cplusplus
}
#endif
