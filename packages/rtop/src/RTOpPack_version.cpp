// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "RTOpPack_version.hpp"
#include "Trilinos_version.h"

std::string RTOpPack::version()
{ 
  return("RTOp in Trilinos " TRILINOS_VERSION_STRING); 
}
