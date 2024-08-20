// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_ABSTRACT_PREPOSTOPERATOR_H
#define NOX_ABSTRACT_PREPOSTOPERATOR_H

#include "NOX_Observer.hpp"

// This file is deprecated but exists for backwards compatibility in
// transition from PrePostOperator to Observer. Please include
// NOX_Observer.hpp instead.

namespace NOX {
  namespace Abstract {
    using PrePostOperator = NOX::Observer;
  }
}

#endif
