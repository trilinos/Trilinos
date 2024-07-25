// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#include "SortingTools.h"


int SortingTools::sortScalars(int n, double *y, int *perm) const {

  // Sort a vector into increasing order of algebraic values
  //
  // Input:
  //
  // n    (integer ) = Size of the array (input)
  // y    (double* ) = Array of length n to be sorted (input/output)
  // perm (integer*) = Array of length n with the permutation (input/output)
  //                   Optional argument

  int i, j;
  int igap = n / 2;

  if (igap == 0) {
    if ((n > 0) && (perm != 0)) {
      perm[0] = 0;
    }
    return 0;
  }

  if (perm) {
    for (i = 0; i < n; ++i)
      perm[i] = i;
  }

  while (igap > 0) {
    for (i=igap; i<n; ++i) {
      for (j=i-igap; j>=0; j-=igap) {
        if (y[j] > y[j+igap]) {
          double tmpD = y[j];
          y[j] = y[j+igap];
          y[j+igap] = tmpD;
          if (perm) {
            int tmpI = perm[j];
            perm[j] = perm[j+igap];
            perm[j+igap] = tmpI;
          }
        }
        else {
          break;
        }
      }
    }
    igap = igap / 2;
  }

  return 0;

}


int SortingTools::sortScalars_Vectors(int num, double *lambda, double *Q, int ldQ) const {

  // This routines sorts the scalars (stored in lambda) in ascending order.
  // The associated vectors (stored in Q) are accordingly ordered.
  // One vector is of length ldQ.
  // Q must be of size ldQ * num.

  int info = 0;
  int i, j;

  int igap = num / 2;

  if ((Q) && (ldQ > 0)) {
    double *vec = new double[ldQ];
    double tmp;
    while (igap > 0) {
      for (i=igap; i < num; ++i) {
        for (j=i-igap; j>=0; j-=igap) {
          if (lambda[j] > lambda[j+igap]) {
            tmp = lambda[j];
            lambda[j] = lambda[j+igap];
            lambda[j+igap] = tmp;
            ////////////////////
            memcpy(vec, Q + j*ldQ, ldQ*sizeof(double));
            memcpy(Q + j*ldQ, Q + (j+igap)*ldQ, ldQ*sizeof(double));
            memcpy(Q + (j+igap)*ldQ, vec, ldQ*sizeof(double));
          } 
          else {
            break;
          }
        }
      }
      igap = igap / 2;
    } // while (igap > 0)
    delete[] vec;
  } // if ((Q) && (ldQ > 0))
  else {
    while (igap > 0) {
      for (i=igap; i < num; ++i) {
        for (j=i-igap; j>=0; j-=igap) {
          if (lambda[j] > lambda[j+igap]) {
            // Swap two scalars
            double tmp = lambda[j];
            lambda[j] = lambda[j+igap];
            lambda[j+igap] = tmp;
          } 
          else {
            break;
          }
        }
      }
      igap = igap / 2;
    } // while (igap > 0)
  } // if ((Q) && (ldQ > 0))

  return info;

}


