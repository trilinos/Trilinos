// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Version.hpp"
#include "Trilinos_version.h"

namespace Tpetra {

  std::string version() {
    return("Tpetra in Trilinos " TRILINOS_VERSION_STRING);
  }

} // namespace Tpetra

