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
  
  - The "Convergence Criteria" is met (the user has two choices): <br>
  1. "Armijo-Goldstein" conditions: \f$ f(x_{n-1}+ \lambda s) \le f(x_{n-1}) +\alpha \lambda f'(x_{n-1}) \f$ <br>
  2. "Ared/Pred" conditions: \f$ f(x_{n-1}+ \lambda s) \le  f(x_{n-1})(1-\alpha(1-\eta)) \f$ where \f$ \eta \f$ is the linear solve tolerance in the inexact Newton method. <br><br>

  - The maximum iterations are reached. This is considered a failure and the recovery step is taken.

  - The minimum step length is reached. This is considered a line search failure and the recovery step is taken.

  The following parameters can be specified for this line search:
 
  - "Default Step" - starting step length (defaults to 1.0)

  - "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)

  - "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")

  - "Max Iters" - maximum number of iterations (i.e., RHS computations)

  - "Convergence Criteria" - choice of algorithm used to determine if the line search algorithm has converged (defaults to "Armijo-Goldstein")

  - "Alpha Factor" - Scaling factor (defaults to 1.0e-4)

  - "Min Bounds Factor" - Factor that limits the minimum size of the new step based on the previous step (defaults to 0.1)

  - "Max Bounds Factor" - Factor that limits the maximum size of the new step based on the previous step (defaults to 0.9)

A sublist for output parameters will be created called "Output" in the parameter list used to instantiate or reset the class.  Valid output parameters are:

  - "Total Number of Line Search Iterations"

  - "Total Number of Failed Line Searches"

  - "Total Number Steps Requiring Line Search" - The total number of Newton steps that could not take a full step and required the line search algorithm to reduce the F-Norm appropriately.
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
  
protected:
  //! sets the output parameters for the linesearch in an "Output" sublist of the "Line Search" sublist. 
  bool setOutputParameters();

private:

  //! Algorithms used to determine convergence of the line search
  enum convergenceCriteriaType {ArmijoGoldstein, AredPred};

  //! Choice of convergence criteria
  convergenceCriteriaType convCriteria;

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
  double minBoundFactor;

  //! Factor that limits the maximum size of the new step based on the previous step.
  double maxBoundFactor;

  //! Pointer to the input parameter list.  We need this to create an "Output" sublist to store output parameters.  
  NOX::Parameter::List* inputList;

  //! Counter for number of times line searching was needed (i.e. number of Newton steps where we could not take a full step length of 1.0)
  int totalNumLSSteps;

  //! Counter for the total number of line search iterations performed
  int totalNumIterations;

  //! Counter for the number of failed line searches
  int totalNumFailedLineSearches;

};
} // namespace LineSearch
} // namespace NOX
#endif
