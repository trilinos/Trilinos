// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TRILINOSCOUPLINGS_VERSION_H
#define TRILINOSCOUPLINGS_VERSION_H

#include "TrilinosCouplings_ConfigDefs.h"
#include "Trilinos_version.h"

string TrilinosCouplings_Version() { 
  return("TrilinosCouplings in Trilinos " TRILINOS_VERSION_STRING); 
}

#endif /* TRILINOSCOUPLINGS_VERSION_H */
