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

#ifndef NOX_LINESEARCH_QUADRATIC_H
#define NOX_LINESEARCH_QUADRATIC_H

#include "NOX_LineSearch_Generic.H" // base class

namespace NOX {
namespace LineSearch {

//! Use a quadratic polynomial fit for the line search.
/*!
  This line search is called via NOX::LineSearch::Manager.
 
  Given the following information at %solver step \f$ n \f$: <br> <br>
  1. Old solution group: \f$ x_{n-1} \f$ <br>
  2. New solution group: \f$ x_n \f$ <br>
  3. Direction vector: \f$ s = x_n - x_{n-1} \f$ <br>

  the quadratic line search attempts to minimize the function:

  \f[ f(x+\lambda s) = \frac{1}{2}||F||_2^2 \f]

  where \f$ \lambda = [0, 1.0] \f$ is the step length.  \f$ \lambda \f$ is computed based upon a quadratic interpolation at the \f$ k^{th} \f$ iterate of the line search by:

  \f[ \lambda_k = \lambda_{k-1} \left[ \frac{-f'(x_{n-1})}{2[f(x_{n-1}+ \lambda_{k-1} s)-f(x_{n-1})-f'(x_{n-1})]} \right] \f]

  The line search iterates over the computation of \f$ \lambda \f$ until:
  
  - The Armijo-Goldstein conditions are met: \f$ f(x_{n-1}+ \lambda s) = f(x_{n-1}) +\alpha \lambda f'(x_{n-1}) \f$

  - The maximum iterations are reached. This is considered a failure and the recovery step is taken.

  - The minimum step length is reached. This is considered a line search failure and the recovery step is taken.

  The following parameters can be specified for this line search:
 
  - "Default Step" - starting step length (defaults to 1.0)

  - "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)

  - "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")

  - "Max Iters" - maximum number of iterations (i.e., RHS computations)

  - "Alpha Factor" - Scaling factor for the Armijo-Goldstein condition (defaults to 1.0e-4)

  - "Bounds Factor" - Factor that limits the minimum size of the new step based on the previous step (defaults to 0.1)

*/

class Quadratic : public Generic {

public:

  //! Constructor
  Quadratic(Parameter::List& params);

  //! Destructor
  ~Quadratic();
  
  // derived
  bool reset(Parameter::List& params);

  // derived
  bool compute(Abstract::Group& newgrp, double& step, 
	       const Abstract::Vector& dir,
	       const Solver::Generic& s);
  
private:

   //! Minimum step length (i.e., when we give up)
  double minStep;
 
 //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  //! Maximum itersions
  int maxIters;

  //! Scaling factor for the Armijo-Goldstein condition
  double alpha;

  //! Factor that limits the minimum size of the new step based on the previous step.
  double boundFactor;

};
} // namespace LineSearch
} // namespace NOX
#endif
