// $Id$
// $Source$

// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_SAFEGUARDED_DIRECTION_HPP
#define NOX_LINESEARCH_SAFEGUARDED_DIRECTION_HPP

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

/* \brief A line search that limits the magnitued of individual entries of the direction update vector.  The limits are defined by a user supplied vector.  NOTE: This is dangerous and not a true line search in the fact that the update direction will be changed on an entry-by-entry basis and may no longer be a descent direction.  Use at your own risk!
*/

class SafeguardedDirection : public Generic {

public:

  //! Constructor
  SafeguardedDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
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
  Teuchos::RCP<NOX::Abstract::Vector> limitDifference_;
  Teuchos::RCP<NOX::Abstract::Vector> scaledUpdate_;

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
