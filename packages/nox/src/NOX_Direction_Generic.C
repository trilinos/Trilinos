// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Direction_Generic.H"
#include "NOX_Solver_LineSearchBased.H"


bool NOX::Direction::Generic::compute(NOX::Abstract::Vector& d,
                      NOX::Abstract::Group& g,
                      const NOX::Solver::LineSearchBased& s)
{
  return compute(d, g, dynamic_cast<const NOX::Solver::Generic&>(s));
}
