// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_LINESEARCH_BACKTRACK_H
#define NOX_LINESEARCH_BACKTRACK_H

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_Abstract_Vector.H" // for NOX::Abstract::Vector::NormType

namespace NOX {
namespace LineSearch {

//! Generic backtracking line search.
/*!  This line search starts with the step length defined by "Default
  Step". It checks to see if the norm of the right hand side (RHS) has
  been reduced. It so, it exits successfully. Otherwise, it reduces
  the step length by a factor of one-half. It continues to repeat this
  procedure until it either finds a reduction in the norm of the RHS
  or the step is less than that specified by "Minimum Step". In the
  later case, the line search has failed, and we take the step defined
  by "Recovery Step".

  This line search can be called via NOX::LineSearch::Manager.
 
  The following parameters can be specified for this line search:
 
<ul>
<li> "Default Step" - starting step length (defaults to 1.0)
<li> "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)
<li> "Recovery Step" - step to take when the line search fails
     (defaults to value for "Default Step")
<li> "Max Iters" - maximum number of iterations (i.e., RHS computations)
<li> "Decrease Condition" - choices are the following:
     <ul>
     <li> "Max Norm" - simple decrease using the infinity norm of the RHS (default)
     <li> "Two Norm" - simple decrease using the Euclidean norm of the RHS
     </ul>
</ul>
 */  

class Backtrack : public Generic {

public:

  //! Constructor
  Backtrack(NOX::Parameter::List& params);

  //! Destructor
  ~Backtrack();
  
  // derived
  bool reset(NOX::Parameter::List& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step, 
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s);
  
private:

  /*! \brief Compute the correct norm (according to the input
      parameters) of the given grp. 
  */
  double getNormF(const NOX::Abstract::Group& grp) const;

private:

   //! Minimum step length (i.e., when we give up)
  double minStep;
 
 //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  //! Maximum number of iterations
  int maxIters;
  
  //! Norm used in simple decrease condition
  NOX::Abstract::Vector::NormType normType;

};
} // namespace LineSearch
} // namespace NOX
#endif
