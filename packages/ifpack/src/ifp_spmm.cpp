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

/*
 * Home Grown Sparse BLAS
 *
 * These are just a subset of the functions described in SPARKER
 * Working Note #3.
 *
 * Would be great if these could be templated some day
 *
 */

#include "Ifpack_config.h"

#include <stdlib.h>
#include <iostream>
using namespace std;
#include "ifp_spblas.h"

#define _SpMatVal(_a,_lda,_row,_col) ((_a)[(_lda)*(_col)+(_row)])

static void CoordMatVec_float(int m, int n, int k, const float &alpha,
        const float *val, const int *indx, const int *jndx,
        const int &nnz,
        const float *b, int ldb, float *c, int ldc)
{
  int i, j;

  // To make the compiler happy
  if (k && m)
    ;

  // Frob these so we can use one-based indexing externally
  b -= 1;
  c -= 1;

  if (alpha == 1.0) {
    if (n == 1)
      for (j = 0; j < nnz; j++)
    c[indx[j]] += b[jndx[j]] * val[j];
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < nnz; j++)
      _SpMatVal(c, ldc, indx[j], i) += _SpMatVal(b, ldb, indx[j], i) * val[j];
  } else {
    if (n == 1)
      for (j = 0; j < nnz; j++)
    c[indx[j]] += alpha * b[jndx[j]] * val[j];
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < nnz; j++)
      _SpMatVal(c, ldc, indx[j], i) +=
        alpha * _SpMatVal(b, ldb, indx[j], i) * val[j];
  }
}

static void CoordMatVec_double(int m, int n, int k, const double &alpha,
        const double *val, const int *indx, const int *jndx,
        const int &nnz,
        const double *b, int ldb, double *c, int ldc)
{
  int i, j;

  // To make the compiler happy
  if (k && m)
    ;

  // Frob these so we can use one-based indexing externally
  b -= 1;
  c -= 1;

  if (alpha == 1.0) {
    if (n == 1)
      for (j = 0; j < nnz; j++)
    c[indx[j]] += b[jndx[j]] * val[j];
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < nnz; j++)
      _SpMatVal(c, ldc, indx[j], i) += _SpMatVal(b, ldb, indx[j], i) * val[j];
  } else {
    if (n == 1)
      for (j = 0; j < nnz; j++)
    c[indx[j]] += alpha * b[jndx[j]] * val[j];
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < nnz; j++)
      _SpMatVal(c, ldc, indx[j], i) +=
        alpha * _SpMatVal(b, ldb, indx[j], i) * val[j];
  }
}

static void
CompColMatVec_double(int m, int n, int k, const double &alpha,
            const double *val, const int *indx, const int *pntr,
            const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;

  if (alpha == 0.0)
    return;

  // To make the compiler happy
  if (m)
    ;

  // Frob these so we can use one-based indexing externally
  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < k; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] += b[i] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < k; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) += _SpMatVal(b, ldb, i, l) * val[j];
  } else {
    if (n == 1)
      for (i = 0; i < k; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] += alpha * b[i] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < k; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) +=
          alpha * _SpMatVal(b, ldb, i, l) * val[j];
  }
}

static void CompColMatVec_float(int m, int n, int k, const float &alpha,
            const float *val, const int *indx, const int *pntr,
            const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;

  if (alpha == 0.0)
    return;

  // To make the compiler happy
  if (m)
    ;

  // Frob these so we can use one-based indexing externally
  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < k; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] += b[i] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < k; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) += _SpMatVal(b, ldb, i, l) * val[j];
  } else {
    if (n == 1)
      for (i = 0; i < k; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] += alpha * b[i] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < k; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) +=
          alpha * _SpMatVal(b, ldb, i, l) * val[j];
  }
}


static void
CompRowMatVec_double(int m, int n, int k, const double &alpha,
            const double *val, const int *indx, const int *pntr,
            const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;

  if (alpha == 0.0)
    return;

  // To make the compiler happy
  if (m || k)
    ;

  // Frob these so we can use one-based indexing externally
  b -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[i] += b[indx[j]] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, i, l) += _SpMatVal(b, ldb, indx[j], l) * val[j];
  } else {
    if (n == 1)
      for (i = 0; i < m; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[i] += alpha * b[indx[j]] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, i, l) +=
          alpha * _SpMatVal(b, ldb, indx[j], l) * val[j];
  }
}


static void
CompRowMatVec_float(int m, int n, int k, const float &alpha,
            const float *val, const int *indx, const int *pntr,
            const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;

  if (alpha == 0.0)
    return;

  // To make the compiler happy
  if (m || k)
    ;

  // Frob these so we can use one-based indexing externally
  b -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[i] += b[indx[j]] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, i, l) += _SpMatVal(b, ldb, indx[j], l) * val[j];
  } else {
    if (n == 1)
      for (i = 0; i < m; i++)
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[i] += alpha * b[indx[j]] * val[j];
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++)
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, i, l) +=
          alpha * _SpMatVal(b, ldb, indx[j], l) * val[j];
  }
}

static void
ScaleRectangularArray_double(int m, int n, double *c, int ldc, 
    const double &beta)
{
  int i, j;

  if (beta == 1.0)
    return;

  if (beta == 0.0) {
    if (n == 1)
      for (j = 0; j < m; j++)
    c[j] = 0.0;
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      _SpMatVal(c, ldc, j, i) = 0.0;
  } else {
    if (n == 1)
      for (j = 0; j < m; j++)
    c[j] *= beta;
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      _SpMatVal(c, ldc, j, i) *= beta;
  }
}

static void
ScaleRectangularArray_float(int m, int n, float *c, int ldc, 
    const double &beta)
{
  int i, j;

  if (beta == 1.0)
    return;

  if (beta == 0.0) {
    if (n == 1)
      for (j = 0; j < m; j++)
    c[j] = 0.0;
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      _SpMatVal(c, ldc, j, i) = 0.0;
  } else {
    if (n == 1)
      for (j = 0; j < m; j++)
    c[j] *= beta;
    else
      for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      _SpMatVal(c, ldc, j, i) *= beta;
  }
}

/*
 * dcoom -- coordinate format matrix-matrix multiply
 *
 * C <- alpha A B + beta C
 *
 * Arguments:
 *
 * int &transa  Indicates how to operate with the sparse matrix
 *      0 : operate with matrix
 *      1 : operate with transpose matrix
 *      2 : operate with conjugate transpose matrix
 *
 * int &m   Number of rows in matrix c
 *
 * int &n   Number of columns in matrix c
 *
 * int &k   Number of rows in matrix b
 *
 * double &alpha Scalar parameter
 *
 * double &beta Scalar parameter
 *
 * int descra[] Descriptor argument.  Nine element integer array
 *      descra[0] matrix structure
 *          0 : general
 *          1 : symmetric
 *          2 : Hermition
 *          3 : Triangular
 *          4 : Anti-Symmetric
 *          5 : Diagonal
 *      descra[1] upper/lower triangular indicator
 *          1 : lower
 *          2 : upper
 *      descra[2] main diagonal type
 *          0 : non-unit
 *          1 : unit
 *      descra[4] repeated indices?
 *          0 : unknown
 *          1 : no repeated indices
 *
 *
 * double *val  scalar array of length nnz containing matrix entries
 *
 * int *indx    integer array of length nnz containing row indices
 *
 * int *jndx    integer array of length nnz containing column indices
 *
 * double *b    rectangular array with first dimension ldb
 *
 * double *c    rectangular array with first dimension ldc
 *
 * double *work scratch array of length lwork.  lwork should be at least
 *      max(m,n)
 *
 */
IFPACK_DEPRECATED void F77NAME(scoomm)
  (const int &transa, const int &m, const int &n, const int &k,
   const float &alpha,
   const int descra[], const float *val,
   const int *indx, const int *jndx, const int &nnz,
   const float *b, const int &ldb,
   const float &beta, float *c, const int &ldc,
   float *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  ScaleRectangularArray_float(m, n, c, ldc, beta);

  if (alpha == 0.0)
    return;

  // Use this hack if transpose is desired
  if (transa == 1 || transa == 2) {
    const int *itmp = indx;
    indx = jndx;
    jndx = itmp;
  }
  CoordMatVec_float(m, n, k, alpha, val, indx, jndx, nnz, b, ldb, c, ldc);
}


IFPACK_DEPRECATED void F77NAME(dcoomm)
  (const int &transa, const int &m, const int &n, const int &k,
   const double &alpha,
   const int descra[], const double *val,
   const int *indx, const int *jndx, const int &nnz,
   const double *b, const int &ldb,
   const double &beta, double *c, const int &ldc,
   double *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  ScaleRectangularArray_double(m, n, c, ldc, beta);

  if (alpha == 0.0)
    return;

  // Use this hack if transpose is desired
  if (transa == 1 || transa == 2) {
    const int *itmp = indx;
    indx = jndx;
    jndx = itmp;
  }
  CoordMatVec_double(m, n, k, alpha, val, indx, jndx, nnz, b, ldb, c, ldc);
}




/*
 * dcscm -- comp sparse column matrix-matrix multiply
 *
 * Arguments:
 *
 * int &transa  Indicates how to operate with the sparse matrix
 *      0 : operate with matrix
 *      1 : operate with transpose matrix
 *      2 : operate with conjugate transpose matrix
 *
 * int &m   Number of rows in matrix c
 *
 * int &n   Number of columns in matrix c
 *
 * int &k   Number of rows in matrix b
 *
 * double &alpha Scalar parameter
 *
 * double &beta Scalar parameter
 *
 * int descra[] Descriptor argument.  Nine element integer array
 *      descra[0] matrix structure
 *          0 : general
 *          1 : symmetric
 *          2 : Hermition
 *          3 : Triangular
 *          4 : Anti-Symmetric
 *          5 : Diagonal
 *      descra[1] upper/lower triangular indicator
 *          1 : lower
 *          2 : upper
 *      descra[2] main diagonal type
 *          0 : non-unit
 *          1 : unit
 *
 * double *val  scalar array of length nnz containing matrix entries
 *
 * int *indx    integer array of length nnz containing row indices
 *
 * int *pntr    integer array of length k+1 such that pntr(j)-pntr(1)
 *      points to location in val of the first element in column j
 *
 * double *b    rectangular array with first dimension ldb
 *
 * double *c    rectangular array with first dimension ldc
 *
 * double *work scratch array of length lwork.  lwork should be at least
 *      max(m,n)
 *
 */
IFPACK_DEPRECATED void F77NAME(scscmm)
  (const int &transa, const int &m, const int &n, const int &k,
   const float &alpha,
   const int descra[], const float *val,
   const int *indx, const int *pntr, const float *b, int &ldb,
   const float &beta, float *c, const int &ldc,
   float *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  ScaleRectangularArray_float(m, n, c, ldc, beta);

  if (transa == 1 || transa == 2)
    CompRowMatVec_float(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
  else
    CompColMatVec_float(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
}


IFPACK_DEPRECATED void F77NAME(dcscmm)
  (const int &transa, const int &m, const int &n, const int &k,
   const double &alpha,
   const int descra[], const double *val,
   const int *indx, const int *pntr, const double *b, int &ldb,
   const double &beta, double *c, const int &ldc,
   double *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  ScaleRectangularArray_double(m, n, c, ldc, beta);

  if (transa == 1 || transa == 2)
    CompRowMatVec_double(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
  else
    CompColMatVec_double(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
}




/*
 * dcsrm -- comp sparse row matrix-matrix multiply
 *
 * Arguments:
 *
 * int &transa  Indicates how to operate with the sparse matrix
 *      0 : operate with matrix
 *      1 : operate with transpose matrix
 *      2 : operate with conjugate transpose matrix
 *
 * int &m   Number of rows in matrix c
 *
 * int &n   Number of columns in matrix c
 *
 * int &k   Number of rows in matrix b
 *
 * double &alpha Scalar parameter
 *
 * double &beta Scalar parameter
 *
 * int descra[] Descriptor argument.  Nine element integer array
 *      descra[0] matrix structure
 *          0 : general
 *          1 : symmetric
 *          2 : Hermition
 *          3 : Triangular
 *          4 : Anti-Symmetric
 *          5 : Diagonal
 *      descra[1] upper/lower triangular indicator
 *          1 : lower
 *          2 : upper
 *      descra[2] main diagonal type
 *          0 : non-unit
 *          1 : unit
 *
 * double *val  scalar array of length nnz containing matrix entries
 *
 * int *indx    integer array of length nnz containing column indices
 *
 * int *pntr    integer array of length k+1 such that pntr(j)-pntr(1)
 *      points to location in val of the first element in row j
 *
 * double *b    rectangular array with first dimension ldb
 *
 * double *c    rectangular array with first dimension ldc
 *
 * double *work scratch array of length lwork.  lwork should be at least
 *      max(m,n)
 *
 */
IFPACK_DEPRECATED void F77NAME(scsrmm)
  (const int &transa, const int &m, const int &n, const int &k,
   const float &alpha,
   const int descra[], const float *val,
   const int *indx, const int *pntr, const float *b, int &ldb,
   const float &beta, float *c, const int &ldc,
   float *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  ScaleRectangularArray_float(m, n, c, ldc, beta);

  if (transa == 1 || transa == 2)
    CompColMatVec_float(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
  else
    CompRowMatVec_float(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
}


IFPACK_DEPRECATED void F77NAME(dcsrmm)
  (const int &transa, const int &m, const int &n, const int &k,
   const double &alpha,
   const int descra[], const double *val,
   const int *indx, const int *pntr, const double *b, int &ldb,
   const double &beta, double *c, const int &ldc,
   double *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  ScaleRectangularArray_double(m, n, c, ldc, beta);

  if (transa == 1 || transa == 2)
    CompColMatVec_double(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
  else
    CompRowMatVec_double(m, n, k, alpha, val, indx, pntr, b, ldb, c, ldc);
}





