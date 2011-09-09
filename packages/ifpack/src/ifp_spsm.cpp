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


/*
 * int &m   Number of rows in matrix c
 *
 * int &n   Number of columns in matrix c
 *
 * int &k   Number of rows in matrix b
 *
 * Assume diagonal elements are in proper place
 *
 * unitd = 1    D = I
 * unitd = 2    left (row scaling)
 * unitd = 3    right (column scaling)
 */

static void
CopyRectangularArray_double(int m, int n,
             const double *b, int ldb, double *c, int ldc)
{
  int i, l;

  if (b == c)
    return;

  if (n == 1)
    for (i = 0; i < m; i++)
      c[i] = b[i];
  else
    for (l = 0; l < n; l++)
      for (i = 0; i < m; i++)
    _SpMatVal(c, ldc, i, l) = _SpMatVal(b, ldb, i, l);
}


static void
CopyRectangularArray_float(int m, int n,
             const float *b, int ldb, float *c, int ldc)
{
  int i, l;

  if (b == c)
    return;

  if (n == 1)
    for (i = 0; i < m; i++)
      c[i] = b[i];
  else
    for (l = 0; l < n; l++)
      for (i = 0; i < m; i++)
    _SpMatVal(c, ldc, i, l) = _SpMatVal(b, ldb, i, l);
}

static void
CompCol_LowerUnitSolve_double(int m, int n, int unitd, const double *dv, 
      double alpha, const double *val, const int *indx, const int *pntr,
      const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_double(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = alpha * c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}


static void
CompCol_LowerUnitSolve_float(int m, int n, int unitd, const float *dv, 
    float alpha, const float *val, const int *indx, const int *pntr,
    const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_float(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = alpha * c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}


static void
CompCol_LowerDiagSolve_double(int m, int n, int unitd, const double *dv, double alpha,
               const double *val, const int *indx, const int *pntr,
               const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_double(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = c[i+pntr[0]] / val[pntr[i]];
    c[i+pntr[0]] = z;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i]];
      _SpMatVal(c, ldc, i+pntr[0], l) = z;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = alpha * c[i+pntr[0]] / val[pntr[i]];
    c[i+pntr[0]] = z;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i]];
      _SpMatVal(c, ldc, i, l) = z;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}

static void
CompCol_LowerDiagSolve_float(int m, int n, int unitd, const float *dv, 
                float alpha,
               const float *val, const int *indx, const int *pntr,
               const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_float(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = c[i+pntr[0]] / val[pntr[i]];
    c[i+pntr[0]] = z;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i]];
      _SpMatVal(c, ldc, i+pntr[0], l) = z;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = alpha * c[i+pntr[0]] / val[pntr[i]];
    c[i+pntr[0]] = z;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i]];
      _SpMatVal(c, ldc, i, l) = z;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}


static void
CompCol_UpperUnitSolve_double(int m, int n, int unitd, const double *dv, 
                double alpha,
               const double *val, const int *indx, const int *pntr,
               const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_double(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];  
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = alpha * c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}


static void
CompCol_UpperUnitSolve_float(int m, int n, int unitd, const float *dv, 
                float alpha,
               const float *val, const int *indx, const int *pntr,
               const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_float(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];  
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = alpha * c[i+pntr[0]];
    for (j = pntr[i]; j < pntr[i+1]; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l);
      for (j = pntr[i]; j < pntr[i+1]; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}


static void
CompCol_UpperDiagSolve_double(int m, int n, int unitd, const double *dv, 
            double alpha,
               const double *val, const int *indx, const int *pntr,
               const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_double(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = c[i+pntr[0]] / val[pntr[i+1]-1];
    c[i+pntr[0]] = z;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i+1]-1];
      _SpMatVal(c, ldc, i+pntr[0], l) = z;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = alpha * c[i+pntr[0]] / val[pntr[i+1]-1];
    c[i+pntr[0]] = z;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i+1]-1];
      _SpMatVal(c, ldc, i+pntr[0], l) = z;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}


static void
CompCol_UpperDiagSolve_float(int m, int n, int unitd, const float *dv, 
                float alpha,
               const float *val, const int *indx, const int *pntr,
               const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  CopyRectangularArray_float(m, n, b, ldb, &c[pntr[0]-1], ldc);

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = c[i+pntr[0]] / val[pntr[i+1]-1];
    c[i+pntr[0]] = z;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i+1]-1];
      _SpMatVal(c, ldc, i+pntr[0], l) = z;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = alpha * c[i+pntr[0]] / val[pntr[i+1]-1];
    c[i+pntr[0]] = z;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      c[indx[j]] -= z * val[j];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = alpha * _SpMatVal(c, ldc, i+pntr[0], l) / val[pntr[i+1]-1];
      _SpMatVal(c, ldc, i+pntr[0], l) = z;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        _SpMatVal(c, ldc, indx[j], l) -= z * val[j];
    }
  }
}


static void
CompRow_LowerUnitSolve_double(int m, int n, int unitd, const double *dv, 
                double alpha,
               const double *val, const int *indx, const int *pntr,
               const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = b[i] - z;
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = (_SpMatVal(b, ldb, i, l) - z);
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+1] = alpha * (b[i] - z);
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = alpha * (_SpMatVal(b, ldb, i, l) - z);
    }
  }
}


static void
CompRow_LowerUnitSolve_float(int m, int n, int unitd, const float *dv, 
                float alpha,
               const float *val, const int *indx, const int *pntr,
               const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = b[i] - z;
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = (_SpMatVal(b, ldb, i, l) - z);
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+1] = alpha * (b[i] - z);
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = alpha * (_SpMatVal(b, ldb, i, l) - z);
    }
  }
}


static void
CompRow_LowerDiagSolve_double(int m, int n, int unitd, const double *dv, 
                double alpha,
               const double *val, const int *indx, const int *pntr,
               const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = (b[i] - z) / val[pntr[i+1]-1];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i+1]-1];
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = alpha * (b[i] - z) / val[pntr[i+1]-1];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        alpha * (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i+1]-1];
    }
  }
}


static void
CompRow_LowerDiagSolve_float(int m, int n, int unitd, const float *dv, 
                float alpha,
               const float *val, const int *indx, const int *pntr,
               const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = (b[i] - z) / val[pntr[i+1]-1];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i+1]-1];
    }
  } else {
    if (n == 1)
      for (i = 0; i < m; i++) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1] - 1; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = alpha * (b[i] - z) / val[pntr[i+1]-1];
      }
    else
      for (l = 0; l < n; l++)
    for (i = 0; i < m; i++) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1] - 1; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        alpha * (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i+1]-1];
    }
  }
}


static void
CompRow_UpperUnitSolve_double(int m, int n, int unitd, const double *dv, double alpha,
               const double *val, const int *indx, const int *pntr,
               const double *b, int ldb, double *c, int ldc)

{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = b[i] - z;
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = (_SpMatVal(b, ldb, i, l) - z);
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = alpha * (b[i] - z);
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = alpha * (_SpMatVal(b, ldb, i, l) - z);
    }
  }
}


static void
CompRow_UpperUnitSolve_float(int m, int n, int unitd, const float *dv, float alpha,
               const float *val, const int *indx, const int *pntr,
               const float *b, int ldb, float *c, int ldc)

{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];

  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = b[i] - z;
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = (_SpMatVal(b, ldb, i, l) - z);
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {
    z = 0;
    for (j = pntr[i]; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = alpha * (b[i] - z);
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {
      z = 0;
      for (j = pntr[i]; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = alpha * (_SpMatVal(b, ldb, i, l) - z);
    }
  }
}


static void
CompRow_UpperDiagSolve_double(int m, int n, int unitd, const double *dv, double alpha,
               const double *val, const int *indx, const int *pntr,
               const double *b, int ldb, double *c, int ldc)
{
  int i, j, l;
  double z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];
  
  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {    
    z = 0;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = (b[i] - z) / val[pntr[i]];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {  
      z = 0;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i]];
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {    
    z = 0;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = alpha * (b[i] - z) / val[pntr[i]];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {  
      z = 0;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        alpha * (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i]];
    }
  }
}


static void
CompRow_UpperDiagSolve_float(int m, int n, int unitd, const float *dv, float alpha,
               const float *val, const int *indx, const int *pntr,
               const float *b, int ldb, float *c, int ldc)
{
  int i, j, l;
  float z;

  // To make the compiler happy
  if (dv)
    ;

  if (unitd != 1) {
    cerr << "unitd != 1 not implemented" << endl;
    exit(1);
  }

  if (alpha == 0.0)
    return;

  c -= 1;
  val -= pntr[0];
  indx -= pntr[0];
  
  if (alpha == 1.0) {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {    
    z = 0;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = (b[i] - z) / val[pntr[i]];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {  
      z = 0;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i]];
    }
  } else {
    if (n == 1)
      for (i = m - 1; i >= 0; i--) {    
    z = 0;
    for (j = pntr[i] + 1; j < pntr[i+1]; j++)
      z += c[indx[j]] * val[j];
    c[i+pntr[0]] = alpha * (b[i] - z) / val[pntr[i]];
      }
    else
      for (l = 0; l < n; l++)
    for (i = m - 1; i >= 0; i--) {  
      z = 0;
      for (j = pntr[i] + 1; j < pntr[i+1]; j++)
        z += _SpMatVal(c, ldc, indx[j], l) * val[j];
      _SpMatVal(c, ldc, i+pntr[0], l) = 
        alpha * (_SpMatVal(b, ldb, i, l) - z) / val[pntr[i]];
    }
  }
}







/*
 * C <- alpha D A^{-1} B + beta C
 * C <- alpha A^{-1} D B + beta C
 */
IFPACK_DEPRECATED void F77NAME(scscsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const float *dv, const float &alpha,
   const int descra[], const float *val,
   const int *indx, const int *pntr, const float *b, int &ldb,
   const float &beta, float *c, const int &ldc,
   float *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  if (beta != 0.0) {
    cerr << "beta != 0 not implemented" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  if (transa == 0) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompCol_LowerDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_LowerUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompCol_UpperDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_UpperUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    }
  } else if (transa == 1 || transa == 2) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompRow_UpperDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_UpperUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompRow_LowerDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_LowerUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else {
      cerr << "Bad descra[1]" << endl;
      exit(1);
    }
  } else {
    cerr << "Bad transa" << endl;
    exit(1);
  }
}


IFPACK_DEPRECATED void F77NAME(scsrsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const float *dv, const float &alpha,
   const int descra[], const float *val,
   const int *indx, const int *pntr, const float *b, int &ldb,
   const float &beta, float *c, const int &ldc,
   float *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  if (beta != 0.0) {
    cerr << "beta != 0 not implemented" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  if (transa == 0) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompRow_LowerDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_LowerUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompRow_UpperDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_UpperUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    }
  } else if (transa == 1 || transa == 2) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompCol_UpperDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_UpperUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompCol_LowerDiagSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_LowerUnitSolve_float(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else {
      cerr << "Bad descra[1]" << endl;
      exit(1);
    }
  } else {
    cerr << "Bad transa" << endl;
    exit(1);
  }
}


IFPACK_DEPRECATED void F77NAME(dcscsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const double *dv, const double &alpha,
   const int descra[], const double *val,
   const int *indx, const int *pntr, const double *b, int &ldb,
   const double &beta, double *c, const int &ldc,
   double *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  if (beta != 0.0) {
    cerr << "beta != 0 not implemented" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  if (transa == 0) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompCol_LowerDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_LowerUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompCol_UpperDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_UpperUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    }
  } else if (transa == 1 || transa == 2) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompRow_UpperDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_UpperUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompRow_LowerDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_LowerUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else {
      cerr << "Bad descra[1]" << endl;
      exit(1);
    }
  } else {
    cerr << "Bad transa" << endl;
    exit(1);
  }
}




IFPACK_DEPRECATED void F77NAME(dcsrsm)
  (const int &transa, const int &m, const int &n,
   const int &unitd, const double *dv, const double &alpha,
   const int descra[], const double *val,
   const int *indx, const int *pntr, const double *b, int &ldb,
   const double &beta, double *c, const int &ldc,
   double *work, const int &lwork)
{
  if (descra[0] != 0) {
    cerr << "Must have general matrix" << endl;
    exit(1);
  }

  if (beta != 0.0) {
    cerr << "beta != 0 not implemented" << endl;
    exit(1);
  }

  // To make the compiler happy
  if (work && lwork)
    ;

  if (transa == 0) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompRow_LowerDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_LowerUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompRow_UpperDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompRow_UpperUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    }
  } else if (transa == 1 || transa == 2) {
    if (descra[1] == 1) {
      if (descra[2] == 0)
    CompCol_UpperDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_UpperUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else if (descra[1] == 2) {
      if (descra[2] == 0)
    CompCol_LowerDiagSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else if (descra[2] == 1)
    CompCol_LowerUnitSolve_double(m, n, unitd, dv, alpha, val, indx, pntr,
                   b, ldb, c, ldc);
      else {
    cerr << "Bad descra[2]" << endl;
    exit(1);
      }
    } else {
      cerr << "Bad descra[1]" << endl;
      exit(1);
    }
  } else {
    cerr << "Bad transa" << endl;
    exit(1);
  }
}




