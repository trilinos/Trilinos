// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_COMMON_H
#define NOX_COMMON_H

#ifndef TRILINOS_NO_CONFIG_H
#include "NOX_Config.h"
#endif

#include "Teuchos_ConfigDefs.hpp"

#define NOX_MAX(a,b) ((a)>(b)) ? (a) : (b)
#define NOX_MIN(a,b) ((a)<(b)) ? (a) : (b)

namespace NOX {
  typedef Teuchos::Ordinal size_type;
}

#endif
