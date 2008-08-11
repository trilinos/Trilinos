/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef _AZ_BLAS_WRAPPERS_H_
#define _AZ_BLAS_WRAPPERS_H_

#include "az_f77func.h"

#if defined(CRAY_T3X)

#define DGEMV_F77   F77_BLAS_MANGLE(sgemv,SGEMV)
#define DGEMM_F77   F77_BLAS_MANGLE(sgemm,SGEMM)
#define DTRSM_F77   F77_BLAS_MANGLE(strsm,STRSM)
#define DTRMM_F77   F77_BLAS_MANGLE(strmm,STRMM)
#define IDAMAX_F77  F77_BLAS_MANGLE(isamax,ISAMAX)
#define DASUM_F77   F77_BLAS_MANGLE(sasum,SASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(saxpy,SAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(scopy,SCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(sdot,SDOT)
#define DSCAL_F77   F77_BLAS_MANGLE(sscal,SSCAL)

#else

#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)

#endif



#ifdef __cplusplus
extern "C" {
#include <stdio.h>
#endif

/* Double precision BLAS 1 */
double PREFIX DASUM_F77(int* n, double x[], int* incx);
void PREFIX DAXPY_F77(int* n, double* alpha, double x[], int* incx, double y[], int* incy);
void PREFIX DCOPY_F77(int* n, double *x, int* incx, double *y, int* incy);
double PREFIX DDOT_F77(int* n, double x[], int* incx, double y[], int* incy);
void PREFIX DSCAL_F77(int* n, double* alpha, double *x, int* incx);
int PREFIX IDAMAX_F77(int* n, double *x, int* incx);


/* Double precision BLAS 2 */
void PREFIX DGEMV_F77(az_fcd, int* m, int* n, double* alpha, double A[], int* lda,
		       double x[], int* incx, double* beta, double y[], int* incy);

/* Double precision BLAS 3 */
void PREFIX DGEMM_F77(az_fcd, az_fcd, int *m, int *
		      n, int *k, double *alpha, double *a, int *lda, 
		      double *b, int *ldb, double *beta, double *c, int *ldc);
void PREFIX DTRMM_F77(az_fcd, az_fcd, az_fcd, az_fcd, 
		      int *m, int *n, double *alpha, double *a, int * lda, double *b, int *ldb);
void PREFIX DTRSM_F77(az_fcd, az_fcd, az_fcd, az_fcd, 
		      int *m, int *n, double *alpha, double *a, int *
		      lda, double *b, int *ldb);

/* Single precision BLAS 1 */
float PREFIX SASUM_F77(int* n, float x[], int* incx);
void PREFIX SAXPY_F77(int* n, float* alpha, float x[], int* incx, float y[], int* incy);
void PREFIX SCOPY_F77(int* n, float *x, int* incx, float *y, int* incy);
float PREFIX SDOT_F77(int* n, float x[], int* incx, float y[], int* incy);
void PREFIX DSCAL_F77(int* n, double* alpha, double *x, int* incx);
int PREFIX IDAMAX_F77(int* n, double *x, int* incx);

/* Single precision BLAS 2 */
void PREFIX SGEMV_F77(az_fcd, int* m, int* n, float* alpha, float A[], int* lda,
		       float x[], int* incx, float* beta, float y[], int* incy);

/* Single precision BLAS 3 */
void PREFIX SGEMM_F77(az_fcd, az_fcd, int *m, int *
		      n, int *k, float *alpha, float *a, int *lda, 
		      float *b, int *ldb, float *beta, float *c, int *ldc);
void PREFIX STRMM_F77(az_fcd, az_fcd, az_fcd, az_fcd, 
		      int *m, int *n, float *alpha, float *a, int * lda, float *b, int *ldb);
void PREFIX STRSM_F77(az_fcd, az_fcd, az_fcd, az_fcd, 
		      int *m, int *n, float *alpha, float *a, int *
		      lda, float *b, int *ldb);


#ifdef __cplusplus
}
#endif

#endif /* _AZ_BLAS_WRAPPERS_H_ */
