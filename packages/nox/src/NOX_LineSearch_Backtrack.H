// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_BACKTRACK_H
#define NOX_LINESEARCH_BACKTRACK_H

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_Abstract_Vector.H" // for NOX::Abstract::Vector::NormType
#include "Teuchos_RCP.hpp"

// Forward declarations
namespace NOX {
  class Utils;
  namespace MeritFunction {
    class Generic;
  }
}

namespace NOX {
namespace LineSearch {

//! Generic backtracking line search.
/*!  This line search starts with the step length defined by "Default
  Step". It checks to see if the norm of the right hand side (RHS) has
  been reduced. If so, it exits successfully. Otherwise, it reduces
  the step length by the reduction factor (defaults to one-half). It
  continues to repeat this procedure until it either finds a reduction
  in the norm of the RHS or the step is less than that specified by
  "Minimum Step". In the later case, the line search has failed, and
  we take the step defined by "Recovery Step".

  This line search can be called via NOX::LineSearch::Manager.

  The following parameters can be specified for this line search in
  the "Backtrack" sublist of the "Line Search" sublist.

<ul>
<li> "Default Step" - starting step length (defaults to 1.0)
<li> "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)
<li> "Recovery Step" - step to take when the line search fails
     (defaults to value for "Default Step")
<li> "Max Iters" - maximum number of iterations (i.e., RHS computations)
<li> "Reduction Factor" - A multiplier between zero and one that reduces the step size between line search iterations
</ul>
 */

class Backtrack : public Generic {

public:

  //! Constructor
  Backtrack(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& params);

  //! Destructor
  ~Backtrack();

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
           const NOX::Abstract::Vector& dir,
           const NOX::Solver::Generic& s);

private:

  //! Printing Utilities
  Teuchos::RCP<NOX::Utils> utils;

  //! Merit function
  Teuchos::RCP<NOX::MeritFunction::Generic> meritFunctionPtr;

   //! Minimum step length (i.e., when we give up)
  double minStep;

 //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  double reductionFactor;

  //! Maximum number of iterations
  int maxIters;

};
} // namespace LineSearch
} // namespace NOX
#endif
