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

#ifndef SORTING_TOOLS_H
#define SORTING_TOOLS_H

#include <cstring>
using std::memcpy;

class SortingTools {

  public:

    int sortScalars(int n, double *y, int *perm = 0) const;

    int sortScalars_Vectors(int,  double *, double * = 0, int = 0) const;

};

#endif
