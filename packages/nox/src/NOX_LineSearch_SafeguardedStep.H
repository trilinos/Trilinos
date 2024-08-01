// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_SAFEGUARDED_STEP_H
#define NOX_LINESEARCH_SAFEGUARDED_STEP_H

#include "NOX_LineSearch_Generic.H" // base class

#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "Teuchos_RCP.hpp"          // class data member

// Forward Declarations
namespace NOX {
  namespace MeritFunction {
    class Generic;
  }
  class LineSearchCounters;
}

namespace NOX {
namespace LineSearch {

/*!
  \brief A line search that determines the step size by capping the magnitue of specific entries of the direction vector.  The limits are specified by a user defined vector.
*/

class SafeguardedStep : public Generic {

public:

  //! Constructor
  SafeguardedStep(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params);

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
           const NOX::Abstract::Vector& dir,
           const NOX::Solver::Generic& s);

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters();

private:

  void printOpeningRemarks() const;

private:

  Teuchos::RCP<NOX::Abstract::Vector> userLimits_;
  Teuchos::RCP<NOX::Abstract::Vector> invLimits_;
  Teuchos::RCP<NOX::Abstract::Vector> scaledUpdate_;
  double lowerStepBound_;
  double upperStepBound_;

  Teuchos::ParameterList* paramsPtr_;
  Teuchos::RCP<Teuchos::ParameterList> validParams_;
  bool useCounter_;
  Teuchos::RCP<NOX::GlobalData> globalDataPtr_;
  NOX::LineSearch::Utils::Printing print_;
  NOX::LineSearchCounters* counter_;

};
} // namespace LineSearch
} // namespace NOX
#endif
