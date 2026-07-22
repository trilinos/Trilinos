// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_Version.hpp"
#include "Trilinos_version.h"

std::string Thyra::Thyra_Version()
{ 
  return("Thyra in Trilinos " TRILINOS_VERSION_STRING); 
}
