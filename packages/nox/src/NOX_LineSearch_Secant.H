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

#ifndef NOX_LINESEARCH_SECANT_H
#define NOX_LINESEARCH_SECANT_H

#include "NOX_LineSearch_Generic.H" // base class

namespace NOX {
namespace LineSearch {

//! Use Secant linesearch adapted from Shewchuk
/*!

  This is a linesearch in which nonlinear iterations for the step length
  (scalar) are performed using a quasi-Newton (secant) method adapted from
  Shewchuk, "An Introduction to the Conjugate Gradient Method Without the
  Agonizing Pain, ed 1 1/4," 1994.

  This linesearch is not currently constrained by conditions such as those of
  Wolfe or Armijo-Goldstein but rather is free to sample any step length in
  hopes of finding a very accurate minimization of the objective function in
  the search direction provided.  The method assumes that locally the negative 
  of the computed residual serves as the gradient of the function to be
  minimized.  Highly accurate linesearch is expected to help convergence
  of nonlinearCG.

  This line search can be called via NOX::LineSearch::Manager.
 
  The following parameters can be specified for this line search:
 
  - "Minimum Step" - minimum step length allowed (when we give up)

  - "Default Step" - starting step length (defaults to 1.0)

  - "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")

  - "Max Iters" - maximum number of iterations (i.e., RHS computations)

  - Initialization of the numerically approximated Hessian is done using
    a backward difference with perturbation (initial step length) of -1.e-5.

  \author Russ Hopper, Org. 9233, Sandia National Labs

*/

class Secant : public Generic {

public:

  //! Constructor
  Secant(Parameter::List& params);

  //! Destructor
  ~Secant();
  
  // derived
  bool reset(Parameter::List& params);

  // derived
  bool compute(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) ;
  
private:

   //! Computed correction to step length; nonlinear iterations for 
   //! step length are assumed to be converged when abs(alpha) < 1.e-8
  double alpha;

   //! Minimum step length (i.e., when we give up)
  double minstep;

   //! Best step length found so far
  double bestStep;
 
 //! Default step
  double defaultstep;

  //! Default step for linesearch failure
  double recoverystep;

  //! Maximum itersions
  int maxiters;

};
} // namespace LineSearch
} // namespace NOX
#endif
