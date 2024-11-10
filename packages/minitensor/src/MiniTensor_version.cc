// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MiniTensor_version.h"
#include "Trilinos_version.h"

std::string minitensor::version()
{ 
  return("MiniTensor in Trilinos " TRILINOS_VERSION_STRING); 
}
