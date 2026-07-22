// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Version.hpp"
#include "Trilinos_version.h"

std::string Piro::Piro_Version()
{ 
  return("Piro in Trilinos " TRILINOS_VERSION_STRING); 
}
