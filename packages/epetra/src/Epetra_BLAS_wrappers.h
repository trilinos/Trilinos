
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

#define PREFIX
#define Epetra_fcd fcd

#define dasum_ SASUM
#define daxpy_ SAXPY
#define dcopy_ SCOPY
#define ddot_ SDOT
#define dnrm2_ SNRM2
#define dscal_ SSCAL
#define idamax_ ISAMAX
#define dgemv_ SGEMV
#define dger_ SGER
#define dtrmv_ STRMV
#define dgemm_ SGEMM
#define dsymm_ SSYMM
#define dtrmm_ STRMM
#define dtrsm_ STRSM

#elif defined(INTEL_CXML)

#define PREFIX __stdcall
#define Epetra_fcd char *, unsigned int

#define sasum_ SASUM
#define saxpy_ SAXPY
#define scopy_ SCOPY
#define sdot_  SDOT
#define snrm2_ SNRM2
#define sscal_ SSCAL
#define isamax_ ISAMAX
#define sgemv_ SGEMV
#define sger_ SGER
#define strmv_ STRMV
#define sgemm_ SGEMM
#define ssymm_ SSYMM
#define strmm_ STRMM
#define strsm_ STRSM

#define dasum_ DASUM
#define daxpy_ DAXPY
#define dcopy_ DCOPY
#define ddot_ DDOT
#define dnrm2_ DNRM2
#define dscal_ DSCAL
#define idamax_ IDAMAX
#define dgemv_ DGEMV
#define dger_ DGER
#define dtrmv_ DTRMV
#define dgemm_ DGEMM
#define dsymm_ DSYMM
#define dtrmm_ DTRMM
#define dtrsm_ DTRSM


#elif defined(INTEL_MKL)

#define PREFIX
#define Epetra_fcd char *

#define sasum_  SASUM
#define saxpy_  SAXPY
#define scopy_  SCOPY
#define sdot_  SDOT
#define snrm2_  SNRM2
#define sscal_  SSCAL
#define isamax_  ISAMAX
#define sgemv_  SGEMV
#define sger_  SGER
#define strmv_  STRMV
#define sgemm_  SGEMM
#define ssymm_  SSYMM
#define strmm_  STRMM
#define strsm_  STRSM

#define dasum_  DASUM
#define daxpy_  DAXPY
#define dcopy_  DCOPY
#define ddot_  DDOT
#define dnrm2_  DNRM2
#define dscal_  DSCAL
#define idamax_  IDAMAX
#define dgemv_  DGEMV
#define dger_  DGER
#define dtrmv_  DTRMV
#define dgemm_  DGEMM
#define dsymm_  DSYMM
#define dtrmm_  DTRMM
#define dtrsm_  DTRSM

#else

/* Define fcd (Fortran Epetra_fcd descriptor) */
#define PREFIX
#define Epetra_fcd char * 

#if defined(__rs6000)
#define dasum_ dasum
#define daxpy_ daxpy
#define dcopy_ dcopy
#define ddot_ ddot
#define dnrm2_ dnrm2
#define dscal_ dscal
#define idamax_ idamax
#define dgemv_ dgemv
#define dger_ dger
#define dtrmv_ dtrmv
#define dgemm_ dgemm
#define dsymm_ dsymm
#define dtrmm_ dtrmm
#define dtrsm_ dtrsm
#endif

#endif

// Double precision BLAS 1
extern "C" double PREFIX dasum_(int* n, double x[], int* incx);
extern "C" void PREFIX daxpy_(int* n, double* alpha, double x[], int* incx, double y[], int* incy);
extern "C" void PREFIX dcopy_(int* n, double *x, int* incx, double *y, int* incy);
extern "C" double PREFIX ddot_(int* n, double x[], int* incx, double y[], int* incy);
extern "C" double PREFIX dnrm2_(int* n, double x[], int* incx);
extern "C" void PREFIX dscal_(int* n, double* alpha, double *x, int* incx);
extern "C" int PREFIX idamax_(int* n, double *x, int* incx);

// Single precision BLAS 1
extern "C" float PREFIX sasum_(int* n, float x[], int* incx);
extern "C" void PREFIX saxpy_(int* n, float* alpha, float x[], int* incx, float y[], int* incy);
extern "C" void PREFIX scopy_(int* n, float *x, int* incx, float *y, int* incy);
extern "C" float PREFIX sdot_(int* n, float x[], int* incx, float y[], int* incy);
extern "C" float PREFIX snrm2_(int* n, float x[], int* incx);
extern "C" void PREFIX sscal_(int* n, float* alpha, float *x, int* incx);
extern "C" int PREFIX isamax_(int* n, float *x, int* incx);

// Double precision BLAS 2
extern "C" void PREFIX dgemv_(Epetra_fcd, int* m, int* n, double* alpha, double A[], int* lda,
		       double x[], int* incx, double* beta, double y[], int* incy);
extern "C" void PREFIX dtrmv_(Epetra_fcd, Epetra_fcd, Epetra_fcd, int *n, 
		      double *a, int *lda, double *x, int *incx);
extern "C" void PREFIX dger_(int *m, int *n, double *alpha, double *x, int *incx, double *y, 
		     int *incy, double *a, int *lda);


// Single precision BLAS 2
extern "C" void PREFIX sgemv_(Epetra_fcd, int* m, int* n, float* alpha, float A[], int* lda,
		       float x[], int* incx, float* beta, float y[], int* incy);
extern "C" void PREFIX strmv_(Epetra_fcd, Epetra_fcd, Epetra_fcd, int *n, 
		      float *a, int *lda, float *x, int *incx);
extern "C" void PREFIX sger_(int *m, int *n, float *alpha, float *x, int *incx, float *y, 
		     int *incy, float *a, int *lda);

// Double precision BLAS 3
extern "C" void PREFIX dgemm_(Epetra_fcd, Epetra_fcd, int *m, int *
		      n, int *k, double *alpha, double *a, int *lda, 
		      double *b, int *ldb, double *beta, double *c, int *ldc);
extern "C" void PREFIX dsymm_(Epetra_fcd, Epetra_fcd, int *m, int * n,
		      double *alpha, double *a, int *lda, 
		      double *b, int *ldb, double *beta, double *c, int *ldc);
extern "C" void PREFIX dtrmm_(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, double *alpha, double *a, int * lda, double *b, int *ldb);
extern "C" void PREFIX dtrsm_(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, double *alpha, double *a, int *
		      lda, double *b, int *ldb);

// Single precision BLAS 3
extern "C" void PREFIX sgemm_(Epetra_fcd, Epetra_fcd, int *m, int *
		      n, int *k, float *alpha, float *a, int *lda, 
		      float *b, int *ldb, float *beta, float *c, int *ldc);
extern "C" void PREFIX ssymm_(Epetra_fcd, Epetra_fcd, int *m, int * n,
		      float *alpha, float *a, int *lda, 
		      float *b, int *ldb, float *beta, float *c, int *ldc);
extern "C" void PREFIX strmm_(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, float *alpha, float *a, int * lda, float *b, int *ldb);
extern "C" void PREFIX strsm_(Epetra_fcd, Epetra_fcd, Epetra_fcd, Epetra_fcd, 
		      int *m, int *n, float *alpha, float *a, int *
		      lda, float *b, int *ldb);

extern "C" void PREFIX xerbla_(Epetra_fcd, int *info);
