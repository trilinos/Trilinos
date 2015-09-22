/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef _AZ_LAPACK_WRAPPERS_H_
#define _AZ_LAPACK_WRAPPERS_H_

#include "az_f77func.h"

#if defined(CRAY_T3X)

#define DGETRF_F77  F77_BLAS_MANGLE(sgetrf,SGETRF)
#define DGETRS_F77  F77_BLAS_MANGLE(sgetrs,SGETRS)
#define DPOTRF_F77  F77_BLAS_MANGLE(spotrf,SPOTRF)
#define DGETRI_F77  F77_BLAS_MANGLE(sgetri,SGETRI)
#define DSTEBZ_F77  F77_BLAS_MANGLE(sstebz,SSTEBZ)
#define DGEEV_F77   F77_BLAS_MANGLE(sgeev,SGEEV)

#else

#define DGETRF_F77  F77_BLAS_MANGLE(dgetrf,DGETRF)
#define DGETRS_F77  F77_BLAS_MANGLE(dgetrs,DGETRS)
#define DPOTRF_F77  F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define DGETRI_F77  F77_BLAS_MANGLE(dgetri,DGETRI)
#define DSTEBZ_F77  F77_BLAS_MANGLE(dstebz,DSTEBZ)
#define DGEEV_F77   F77_BLAS_MANGLE(dgeev,DGEEV)

#endif

#ifdef __cplusplus
extern "C" {
#include <stdio.h>
#endif


  /* Double precision LAPACK linear solvers */
void PREFIX DGETRF_F77(int* m, int* n, double* a, int* lda, int* ipiv, int* info); 
void PREFIX DGETRS_F77(az_fcd, int* n, int* nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void PREFIX DGETRI_F77(int* n, double* a, int* lda, int*ipiv, double * work , int* lwork, int* info);
void PREFIX DPOTRF_F77(az_fcd, int* n, double* a, int* lda, int* info);
void PREFIX DSTEBZ_F77(az_fcd, az_fcd, int *, double *, double *, int *, int *,
		       double *, double *, double *, int *, int *, double *, int *,
		       int *, double *, int *, int *);
void PREFIX DGEEV_F77(az_fcd, az_fcd, int *, double *, int *, double *,
		      double *, double *, int *, double *, int *,
		      double *, int *, int *);

  /* Single precision LAPACK linear solvers*/
void PREFIX SGETRF_F77(int* m, int* n, float* a, int* lda, int* ipiv, int* info); 
void PREFIX SGETRS_F77(az_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
void PREFIX SGETRI_F77(int* n, float* a, int* lda, int*ipiv, float * work , int* lwork, int* info);
void PREFIX SPOTRF_F77(az_fcd, int* n, float* a, int* lda, int* info); 
void PREFIX SSTEBZ_F77(az_fcd, az_fcd, int *, float *, float *, int *, int *,
		       float *, float *, float *, int *, int *, double*,  int *,
		       int *, float *, int *, int *);
void PREFIX SGEEV_F77(az_fcd, az_fcd, int *, float *, int *, float *,
		      float *, float *, int *, float *, int *,
		      float *, int *, int *);


#ifdef __cplusplus
}
#endif

#endif /* _AZ_LAPACK_WRAPPERS_H_ */

void PREFIX DGETRS_F77(az_fcd, int* n, int* nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void PREFIX SGETRS_F77(az_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
