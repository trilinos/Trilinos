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

#ifndef _AZ_LAPACK_WRAPPERS_H_
#define _AZ_LAPACK_WRAPPERS_H_

#include "az_f77func.h"

#if defined(CRAY_T3X)

#define DGETRF_F77  F77_FUNC(sgetrf,SGETRF)
#define DGETRS_F77  F77_FUNC(sgetrs,SGETRS)
#define DPOTRF_F77  F77_FUNC(spotrf,SPOTRF)
#define DGETRI_F77  F77_FUNC(sgetri,SGETRI)
#define DLASWP_F77  F77_FUNC(slaswp,SLASWP)
#define DLAIC1_F77  F77_FUNC(slaic1,SLAIC1)

#else

#define DGETRF_F77  F77_FUNC(dgetrf,DGETRF)
#define DGETRS_F77  F77_FUNC(dgetrs,DGETRS)
#define DPOTRF_F77  F77_FUNC(dpotrf,DPOTRF)
#define DGETRI_F77  F77_FUNC(dgetri,DGETRI)
#define DLASWP_F77  F77_FUNC(dlaswp,DLASWP)
#define DLAIC1_F77  F77_FUNC(dlaic1,DLAIC1)

#endif

#ifdef __cplusplus
extern "C" {
#include <stdio.h>
#endif

void DLASWP_F77(int *, double *, int *, int *, int *, int *, int *);

void DLAIC1_F77(int * , int *, double *, double *, double *, double *,
			  double *, double *, double *);
void SLASWP_F77(int *, float *, int *, int *, int *, int *, int *);

void SLAIC1_F77(int * , int *, float *, float *, float *, float *,
			  float *, float *, float *);

  /* Double precision LAPACK linear solvers */
void DGETRF_F77(int* m, int* n, double* a, int* lda, int* ipiv, int* info); 
void DGETRS_F77(az_fcd, int* n, int* nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void DGETRI_F77(int* n, double* a, int* lda, int*ipiv, double * work , int* lwork, int* info);
void DPOTRF_F77(az_fcd, int* n, double* a, int* lda, int* info); 

  /* Single precision LAPACK linear solvers*/
void SGETRF_F77(int* m, int* n, float* a, int* lda, int* ipiv, int* info); 
void SGETRS_F77(az_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
void SGETRI_F77(int* n, float* a, int* lda, int*ipiv, float * work , int* lwork, int* info);
void SPOTRF_F77(az_fcd, int* n, float* a, int* lda, int* info); 


#ifdef __cplusplus
}
#endif

#endif /* _AZ_LAPACK_WRAPPERS_H_ */

void DGETRS_F77(az_fcd, int* n, int* nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void SGETRS_F77(az_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
