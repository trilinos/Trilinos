// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_PREPOSTOPERATOR_VECTOR_H
#define NOX_PREPOSTOPERATOR_VECTOR_H

#include "NOX_Observer_Vector.hpp"

// This file is deprecated, but exists for backwards compatibility to
// support transition from PrePostOperator to Observer. Please include
// NOX_Observer_Vector.hpp instead.

namespace NOX {
  using PrePostOperatorVector = NOX::ObserverVector;
}

#endif
