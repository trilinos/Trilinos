// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_LINESEARCH_BACKTRACK_H
#define NOX_LINESEARCH_BACKTRACK_H

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_Abstract_Vector.H" // for NOX::Abstract::Vector::NormType
#include "Teuchos_RefCountPtr.hpp"

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
  Backtrack(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	    Teuchos::ParameterList& params);

  //! Destructor
  ~Backtrack();
  
  // derived
  bool reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	     Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step, 
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s);
  
private:
  
  //! Printing Utilities
  Teuchos::RefCountPtr<NOX::Utils> utils;

  //! Merit function
  Teuchos::RefCountPtr<NOX::MeritFunction::Generic> meritFunctionPtr;

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
