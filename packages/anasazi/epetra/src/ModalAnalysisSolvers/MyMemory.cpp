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

#include "MyMemory.h"

char *startingPoint;

void initMemCounters() {

#ifdef INTEL_CXML
  startingPoint = NULL;
#else
  startingPoint = (char *) sbrk(0);
#endif

  return;

}


double currentSize() {

#ifdef INTEL_CXML
  return 0.0;
#else
  char *current = (char *) sbrk(0);
  return (current - startingPoint)/(1024.0*1024.0);
#endif

}

