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

#ifndef NOX_LINESEARCH_POLYNOMIAL_H
#define NOX_LINESEARCH_POLYNOMIAL_H

#include "NOX_LineSearch_Generic.H" // base class

namespace NOX {
namespace LineSearch {

//! Use quad/cubic polynomial fit from Dennis & Schnabel Sec. 6.3.2.
/*!
  This line search is adapted from "Numerical Methods for Unconstrained
  Optimization and Nonlinear Equations," J. E. Dennis Jr. and Robert B.  
  Schnabel, Sec. 6.3.2 and can be called via NOX::LineSearch::Manager.
 
  Given the following information at %solver step \f$ n \f$: <br> <br>
  1. Old solution group: \f$ x_{n-1} \f$ <br>
  2. New solution group: \f$ x_n \f$ <br>
  3. Direction vector: \f$ s = x_n - x_{n-1} \f$ <br>

  the linesearch (backtracking) attempts to minimize the function :

  \f[ f(x+\lambda s) = \frac{1}{2}||F||_2^2 \f]

  where \f$ \lambda = [0, 1.0] \f$ is the step length.  \f$ \lambda
  \f$ is computed based on polynomial fits to information on \f$ f \f$
  and its derivative \f$ f' \f$ at the old solution and newly computed
  values of \f$ f \f$.

  The polynomial linesearch starts by setting \f$ \lambda_0 \f$
  to the value specified by "Default Step" (see below) and checking
  acceptability criteria (see below).  If the step is unacceptable,
  a new step is computed as the minimum to a quadratic fit resulting in:

  \f[ \lambda_1 = \frac{-f'(x_{n-1})}{2\left[f(x_{n-1}+\lambda_0 s)-f(x_{n-1})
      -f'(x_{n-1})\right]} \f]

  where \f$ f'(x_{n-1}) \equiv \left. J^TF\right|_{x_{n-1}} \f$.

  If this step is still unacceptable, new steps are computed as minimizations
  to cubic fits with the following result for \f$ k>1 \f$ :

  \f[ \lambda_k = \frac{-b+\sqrt{b^2-3af'(x_{n-1})}}{3a} \f]

  where

  \f[ \left[ \begin{array}{c} a \\ \\ b \end{array} \right] = 
      \frac{1}{\lambda_{k-1} - \lambda_{k-2}} \left[ \begin{array}{cc}
      \lambda_{k-1}^{-2} & -\lambda_{k-2}^{-2} \\ & \\
      -\lambda_{k-2}\lambda_{k-1}^{-2} & \lambda_{k-1}\lambda_{k-2}^{-2}
      \end{array} \right]
      \left[ \begin{array}{c} f(x_{n-1}+\lambda_{k-1} s) - f(x_{n-1}) -
                              f'(x_{n-1})\lambda_{k-1} \\ \\
                              f(x_{n-1}+\lambda_{k-2} s) - f(x_{n-1}) -
                              f'(x_{n-1})\lambda_{k-2} \end{array} \right]
  \f]

  <B>Step Length Termination </B> <br>
  The line search iterates over the computation of \f$ \lambda \f$ until:

  - The "Convergence Criteria" is met (the user has two choices): <br>
  1. "Armijo-Goldstein" conditions: \f$ f(x_{n-1}+ \lambda s) \le f(x_{n-1}) +\alpha \lambda f'(x_{n-1}) \f$ <br>
  2. "Ared/Pred" conditions: \f$ f(x_{n-1}+ \lambda s) \le  f(x_{n-1})(1-\alpha(1-\eta)) \f$ where \f$ \eta \f$ is the linear solve tolerance in the inexact Newton method. <br><br>

  - The maximum iterations are reached. This is considered a failure
    and the recovery step is taken.

  - The minimum step length is reached. This is considered a line search
    failure and the recovery step is taken.

  <B> Step bounds </B> <br> 
  The following bounds are currently enforced on the step length calculation:

  - For \f$ k>0 \f$, \f$ \lambda_k \ge \gamma_{min} \lambda_{k-1} \f$

  - For \f$ k>1 \f$, \f$ \lambda_k \le \gamma_{max} \lambda_{k-1} \f$

  where \f$ \gamma_{min} \f$ and \f$ \gamma_{max} \f$ are bounds factors that limit the change in step size between the iterations of the step length calculation (\f$ k \f$). 

  <B> Input Parameters </B> <br>
  The following parameters can be specified for this line search:
 
  - "Default Step" - starting step length (defaults to 1.0)

  - "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)

  - "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")

  - "Max Iters" - maximum number of line search iterations (default to 100)

  - "Convergence Criteria" - choice of algorithm used to determine if the line search algorithm has converged (defaults to "Armijo-Goldstein").  Choices are (See above for explanations): <br>
  1. "Armijo-Goldstein" <br>
  2. "Ared/Pred" <br><br>

  - "Alpha Factor" - Scaling factor, \f$ \alpha \f$ (defaults to 1.0e-4)

  - "Min Bounds Factor" - \f$ \gamma_{min} \f$, the factor that limits the minimum size of the new step based on the previous step (defaults to 0.1)

  - "Max Bounds Factor" - \f$ \gamma_{max} \f$, the factor that limits the maximum size of the new step based on the previous step (defaults to 0.9)

  <B> Output Parameters </B> <br>
A sublist for output parameters will be created called "Output" in the parameter list used to instantiate or reset the class.  Valid output parameters are:

  - "Total Number of Line Search Calls" - Total number of calls to the compute() method of this line search.

  - "Total Number of Non-trivial Line Searches" - The total number of steps that could not directly take a full step and meet the required "Convergence Criteria" (i.e. The line search had to reduce the step length using inner iteration calculations over iterate \f$ k \f$). 

  - "Total Number of Failed Line Searches" - total number of line searches that failed and used a recovery step.

  - "Total Number of Line Search Inner Iterations" - total number of inner iterations \f$ k \f$ performed by this object.

  \author Russ Hooper, Org. 9233, Sandia National Labs
*/

class Polynomial : public Generic {

public:

  //! Constructor
  Polynomial(Parameter::List& params);

  //! Destructor
  ~Polynomial();
  
  // derived
  bool reset(Parameter::List& params);

  // derived
  bool compute(Abstract::Group& newgrp, double& step, 
	       const Abstract::Vector& dir,
	       const Solver::Generic& s);
	       
protected:

  //! sets the output parameters for the linesearch in an "Output" sublist of the "Line Search" sublist. 
  bool setOutputParameters();  
  
protected:

  //! Algorithms used to determine convergence of the line search
  enum convergenceCriteriaType {ArmijoGoldstein, AredPred};

  //! Choice of convergence criteria
  convergenceCriteriaType convCriteria;

   //! Minimum step length (i.e., when we give up)
  double minstep;
 
 //! Default step
  double defaultstep;

  //! Default step for linesearch failure
  double recoverystep;

  //! Maximum itersions
  int maxiters;

  //! Scaling factor for the Armijo-Goldstein condition
  double alpha;

  //! Factor that limits the minimum size of the new step based on the previous step.
  double minBoundFactor;

  //! Factor that limits the maximum size of the new step based on the previous step.
  double maxBoundFactor;

  //! Pointer to the input parameter list.  We need this to create an "Output" sublist to store output parameters.  
  NOX::Parameter::List* inputList;

  //! Counter for the total number of times the compute() method was called on this object. 
  int totalNumLineSearchCalls;

  //! Counter for number of times line searching was needed (i.e. number of Newton steps where we could not take a full step length).
  int totalNumNonTrivialLineSearches;

  //! Counter for the number of failed line searches.
  int totalNumFailedLineSearches;

  //! Counter for the total number of inner iterations performed.
  int totalNumIterations;

};
} // namespace LineSearch
} // namespace NOX
#endif
