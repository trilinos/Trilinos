// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_FULLSTEP_H
#define NOX_LINESEARCH_FULLSTEP_H

#include "NOX_LineSearch_Generic.H" // base class

namespace NOX {
namespace LineSearch {

//! Simplest line search - always take the full step.
/*!
  This line search can be called via NOX::LineSearch::Manager.

  The following parameters can be specified in the "Full Step" sublist
  of the "Line Search" sublist:

  - "Full Step" - length of a full step (defaults to 1.0)
*/

class FullStep : public Generic {

public:

  //! Constructor
  FullStep(const Teuchos::RCP<NOX::GlobalData>& gd,
       Teuchos::ParameterList& params);

  //! Destructor
  ~FullStep();

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
           const NOX::Abstract::Vector& dir,
           const NOX::Solver::Generic& s);
private:

  //! Full step length
  double stepSize;
};
} // namespace LineSearch
} // namespace NOX

#endif
