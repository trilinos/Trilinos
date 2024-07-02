// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_VERSION_CPP
#define AMESOS2_VERSION_CPP

#include <string>

#include "Amesos2_config.h"
#include "Trilinos_version.h"

namespace Amesos2 {

std::string version() {
  return("Amesos2 version " AMESOS2_VERSION " in Trilinos " TRILINOS_VERSION_STRING);
}

} // namespace Amesos2

#endif
