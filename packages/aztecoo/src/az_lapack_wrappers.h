/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef _AZ_LAPACK_WRAPPERS_H_
#define _AZ_LAPACK_WRAPPERS_H_

#include "az_f77func.h"

#if defined(CRAY_T3X)

#define DGETRF_F77  F77_FUNC(sgetrf,SGETRF)
#define DGETRS_F77  F77_FUNC(sgetrs,SGETRS)
#define DPOTRF_F77  F77_FUNC(spotrf,SPOTRF)
#define DGETRI_F77  F77_FUNC(sgetri,SGETRI)

#else

#define DGETRF_F77  F77_FUNC(dgetrf,DGETRF)
#define DGETRS_F77  F77_FUNC(dgetrs,DGETRS)
#define DPOTRF_F77  F77_FUNC(dpotrf,DPOTRF)
#define DGETRI_F77  F77_FUNC(dgetri,DGETRI)

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

  /* Single precision LAPACK linear solvers*/
void PREFIX SGETRF_F77(int* m, int* n, float* a, int* lda, int* ipiv, int* info); 
void PREFIX SGETRS_F77(az_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
void PREFIX SGETRI_F77(int* n, float* a, int* lda, int*ipiv, float * work , int* lwork, int* info);
void PREFIX SPOTRF_F77(az_fcd, int* n, float* a, int* lda, int* info); 


#ifdef __cplusplus
}
#endif

#endif /* _AZ_LAPACK_WRAPPERS_H_ */

void PREFIX DGETRS_F77(az_fcd, int* n, int* nrhs, double* a,
                       int* lda, int*ipiv, double*x , int* ldx, int* info);
void PREFIX SGETRS_F77(az_fcd, int* m, int* n, float* a,
                       int* lda, int*ipiv, float*x , int* ldx, int* info);
