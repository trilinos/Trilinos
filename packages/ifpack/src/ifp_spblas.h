/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5         */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef IFP_SPBLAS_H
#define IFP_SPBLAS_H

#include "Ifpack_config.h"

#include "ifp_arch.h"

#ifdef COMPLEX_SUPPORT
#include "complex.h"
#endif

extern "C" {

IFPACK_DEPRECATED void F77NAME(scoomm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const float &alpha, 
   const int descra[], const float *val, 
   const int *indx, const int *jndx, const int &nnz, 
   const float *b, const int &ldb, 
   const float &beta, float *c, const int &ldc, 
   float *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(scscmm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const float &alpha, 
   const int descra[], const float *val, 
   const int *indx, const int *pntr, const float *b, int &ldb, 
   const float &beta, float *c, const int &ldc, 
   float *work, const int &lwork);
   
IFPACK_DEPRECATED void F77NAME(scsrmm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const float &alpha, 
   const int descra[], const float *val, 
   const int *indx, const int *pntr, const float *b, int &ldb, 
   const float &beta, float *c, const int &ldc, 
   float *work, const int &lwork);
   
IFPACK_DEPRECATED void F77NAME(dcoomm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const double &alpha, 
   const int descra[], const double *val, 
   const int *indx, const int *jndx, const int &nnz, 
   const double *b, const int &ldb, 
   const double &beta, double *c, const int &ldc, 
   double *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(dcscmm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const double &alpha, 
   const int descra[], const double *val, 
   const int *indx, const int *pntr, const double *b, int &ldb, 
   const double &beta, double *c, const int &ldc, 
   double *work, const int &lwork);
   
IFPACK_DEPRECATED void F77NAME(dcsrmm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const double &alpha, 
   const int descra[], const double *val, 
   const int *indx, const int *pntr, const double *b, int &ldb, 
   const double &beta, double *c, const int &ldc, 
   double *work, const int &lwork);
  

IFPACK_DEPRECATED void F77NAME(dcscsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const double *dv, const double &alpha, 
   const int descra[], const double *val, 
   const int *indx, const int *pntr, const double *b, int &ldb, 
   const double &beta, double *c, const int &ldc, 
   double *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(dcsrsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const double *dv, const double &alpha, 
   const int descra[], const double *val, 
   const int *indx, const int *pntr, const double *b, int &ldb, 
   const double &beta, double *c, const int &ldc, 
   double *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(scscsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const float *dv, const float &alpha, 
   const int descra[], const float *val, 
   const int *indx, const int *pntr, const float *b, int &ldb, 
   const float &beta, float *c, const int &ldc, 
   float *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(scsrsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const float *dv, const float &alpha, 
   const int descra[], const float *val, 
   const int *indx, const int *pntr, const float *b, int &ldb, 
   const float &beta, float *c, const int &ldc, 
   float *work, const int &lwork);

#ifdef COMPLEX_SUPPORT

IFPACK_DEPRECATED void F77NAME(zcoomm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const complex &alpha, 
   const int descra[], const complex *val, 
   const int *indx, const int *jndx, const int &nnz, 
   const complex *b, const int &ldb, 
   const complex &beta, complex *c, const int &ldc, 
   complex *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(zcscmm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const complex &alpha, 
   const int descra[], const complex *val, 
   const int *indx, const int *pntr, const complex *b, int &ldb, 
   const complex &beta, complex *c, const int &ldc, 
   complex *work, const int &lwork);
   
IFPACK_DEPRECATED void F77NAME(zcsrmm)
  (const int &transa, const int &m, const int &n, const int &k, 
   const complex &alpha, 
   const int descra[], const complex *val, 
   const int *indx, const int *pntr, const complex *b, int &ldb, 
   const complex &beta, complex *c, const int &ldc, 
   complex *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(zcscsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const complex *dv, const complex &alpha, 
   const int descra[], const complex *val, 
   const int *indx, const int *pntr, const complex *b, int &ldb, 
   const complex &beta, complex *c, const int &ldc, 
   complex *work, const int &lwork);

IFPACK_DEPRECATED void F77NAME(zcsrsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const complex *dv, const complex &alpha, 
   const int descra[], const complex *val, 
   const int *indx, const int *pntr, const complex *b, int &ldb, 
   const complex &beta, complex *c, const int &ldc, 
   complex *work, const int &lwork);


#endif
// COMPLEX_SUPPORT



}

#endif
