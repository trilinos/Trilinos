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

#include "NOX_LineSearch_Common.H" // base class

namespace NOX {
namespace LineSearch {

//! Use quad/cubic polynomial fit from Dennis & Schnabel Sec. 6.3.2.
/*!
  
  This line search is based on materials in the following books:

  - Section 6.3.2 and Algorithm 6.3.1 of "Numerical Methods for
    Unconstrained Optimization and Nonlinear Equations," J. E. Dennis
    Jr. and Robert B.  Schnabel, Prentice Hall, 1983.

  - Section 3.4 of "Numerical Optimization," Jorge Nocedal and Stephen
    J. Wright, Springer, 1999.

  The line search can be called via NOX::LineSearch::Manager.
 
  We are given the following information at %solver step \f$ n \f$: 
  -# Old solution group: \f$ x_{\rm old} \f$ 
  -# Direction vector: \f$ s \f$ 

  The line search attempts to minimize the function 
  \f[ 
  \phi(\lambda) \equiv \frac{1}{2}||F (x_{\rm old} + \lambda s)||_2^2 
  \f]
  where \f$ \lambda \f$ is the step length.  

  We compute \f$ \lambda \f$ based on polynomial fits to information
  on \f$ \phi \f$ and its derivative \f$ \phi' \f$ at the old solution and
  newly computed values of \f$ \phi \f$.

  The polynomial linesearch starts by setting \f$ \lambda_0 \f$
  to the value specified by "Default Step" (see below).

  We generally check that the default step satisfies the acceptability
  criteria (see below) before doing any interpolation.  However, that
  behavior can be overridden by setting "Force Interpolation" to
  true.
  
  The first interpolation step is always quadratic.  Depending on the
  choice of "Interpolation Type" (see below), the following steps are
  either quadratic or cubic.
  
  If we are using a \b quadratic \b model, the \f$k\f$-th step approximation is
  calculated as

  \f[ 
  \lambda_{k} = \frac{-\phi'(0) \lambda_{k-1}^2 }{2 \left[ \phi(\lambda_{k-1}) - \phi(0) 
  -\phi'(0) \lambda_{k-1} \right]} 
  \f]

  where \f$ \phi'(0) = s^T \left(J_{\rm old}^T F_{\rm old} \right) \f$.

  If we are using a \b cubic \b model, the \f$k\f$-th step (\f$ k>1
  \f$) approximation is calculated as

  \f[ 
  \lambda_k = \frac{-b+\sqrt{b^2-3a \phi'(0)}}{3a} 
  \f]

  where

  \f[
      \left[ \begin{array}{c} a \\ \\ b \end{array} \right] = 
      \frac{1}{\lambda_{k-1} - \lambda_{k-2}} \left[ \begin{array}{cc}
      \lambda_{k-1}^{-2} & -\lambda_{k-2}^{-2} \\ & \\
      -\lambda_{k-2}\lambda_{k-1}^{-2} & \lambda_{k-1}\lambda_{k-2}^{-2}
      \end{array} \right]
      \left[ \begin{array}{c} \phi(\lambda_{k-1}) - \phi(0) -
                              \phi'(0)\lambda_{k-1} \\ \\
                              \phi(\lambda_{k-2}) - \phi(0) -
                              \phi'(0)\lambda_{k-2} \end{array} \right]
  \f]

  <b>Step Length Termination</b>

  The line search iterates over the computation of \f$ \lambda \f$ until:

- The "Convergence Criteria" is met (the user has two choices): 
  -# "Armijo-Goldstein" conditions: \f$ f(x_{n-1}+ \lambda s) \le
     f(x_{n-1}) +\alpha \lambda \phi'(0) \f$ 
  -# "Ared/Pred" conditions: \f$ f(x_{n-1}+ \lambda s) \le
     f(x_{n-1})(1-\alpha(1-\eta)) \f$ where \f$ \eta \f$ is the linear
     solve tolerance in the inexact Newton method. 
  -# "None" - If "Force Interpolation" is true, accept the result of
     the first interpolation. Otherwise, accept the default step.
- The maximum iterations are reached. This is considered a failure
  and the recovery step is taken.
- The minimum step length is reached. This is considered a line search
  failure and the recovery step is taken.

  <B> Step bounds </B> 
 
  The following bounds are currently enforced on the step length calculation:

  - For \f$ k>0 \f$, \f$ \lambda_k \ge \gamma_{min} \lambda_{k-1} \f$

  - For \f$ k>1 \f$, \f$ \lambda_k \le \gamma_{max} \lambda_{k-1} \f$ 

  where \f$ \gamma_{min} \f$ and \f$ \gamma_{max} \f$ are bounds
  factors that limit the change in step size between the iterations of
  the step length calculation (\f$ k \f$).

  <B> Input Parameters </B> 

  The following parameters can be specified in the "Polynomial"
  sublist of the "Line Search" sublist.
 
  - "Default Step" - starting step length (defaults to 1.0)

  - "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)

  - "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")

  - "Max Iters" - maximum number of line search iterations (default to 100)

  - "Convergence Criteria" - choice of algorithm used to determine if
    the line search algorithm has converged (defaults to
    "Armijo-Goldstein").  Choices are (See above for explanations):

    -# "Armijo-Goldstein" 
    -# "Ared/Pred" 
    -# "None"

  - "Alpha Factor" - Used in Armijo-Goldstein stopping criteria \f$
    \alpha \f$ (defaults to 1.0e-4)

  - "Min Bounds Factor" - \f$ \gamma_{min} \f$, the factor that limits
    the minimum size of the new step based on the previous step
    (defaults to 0.1)

  - "Max Bounds Factor" - \f$ \gamma_{max} \f$, the factor that limits
    the maximum size of the new step based on the previous step
    (defaults to 0.9)

  - "Force Interpolation" (boolean) - Set to true at least one
    interpolation step should be used. The default is false which
    means that the line search will stop if the default step length
    satisfies the convergence criteria. (defaults to false)

  - "Interpolation Type" - Type of interpolation that should be used. Choices are
    
     <ol>
     <li> "Cubic" (Quadratic for the first iteration)
     <li> "Quadratic"
     </ol>

  <B> Output Parameters </B> 

A sublist for output parameters will be created called "Output" in the
parameter list used to instantiate or reset the class.  Valid output
parameters are:

  - "Total Number of Line Search Calls" - Total number of calls to the compute() method of this line search.

  - "Total Number of Non-trivial Line Searches" - The total number of
    steps that could not directly take a full step and meet the
    required "Convergence Criteria" (i.e. The line search had to
    reduce the step length using inner iteration calculations over
    iterate \f$ k \f$).

  - "Total Number of Failed Line Searches" - total number of line
    searches that failed and used a recovery step.

  - "Total Number of Line Search Inner Iterations" - total number of
    inner iterations \f$ k \f$ performed by this object.

  <B> Flexible step acceptance capability </B> <br>

Line searches typically allow steps to be accepted only if they meet sufficient decrease and curvature criteria, eg Armijo-Goldstein conditions.  It can sometimes happen that such a requirement leads to very many nonlinear iterations being taken before the solution falls into a much steeper basin of attraction.  The following two parameters allow the user to relax the curvature and decrease conditions.  !! WARNING !! Convergence proofs and other guarantees associated with this Line Search no longer pertain when using these options.  

  - "Allowed Relative Increase" - (default 100) - the allowed value \f$ k > 0 \f$ for which 
  \f$
  \frac{||F (x_{\rm old} + \lambda s)||_2}{||F (x_{\rm old})||_2} \le k
  \f$

  - "Maximum Increase Steps" - (defaults to "Max Iters") - maximum number of outer iterations for which the preceding acceptability criterion is allowed.  After this number of outer iterations, the original curvature and decrease conditions must be met for a step to be accepted.

  \author Russ Hooper, Org. 9233, Sandia National Labs
*/

class Polynomial : public Common {

public:

  //! Constructor
  Polynomial(const NOX::Utils& u, NOX::Parameter::List& params);

  //! Destructor
  ~Polynomial();
  
  // derived
  bool reset(NOX::Parameter::List& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step, 
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s);
	       
protected:

  //! Check for convergence
  /*!
    In the case of "Armijo-Goldstein" criteria, we return true if
    \f[
    f_{\rm new} \leq f_{\rm old} + \alpha \; \lambda \; {\rm slope}.
    \f]
    Here \f$\lambda\f$ represents the step and the slope should be 
    \f$d^T \nabla f(x_{\rm old})\f$ where \f$d\f$ is the search direction.
    The slope can be computed via NOX::LineSearch::Common::computeSlope.

    In the case of "Ared/Pred", we return true if
    \f[
    f_{\rm new} \leq f_{\rm old} + (1 - \alpha (1 - \eta))
    \f]
  */
  bool isSufficientDecrease(double newf, double oldf, double step, double slope, double eta) const;

  //! Provides a check to see if a step leading to an increase in the residual is allowed
  bool isIncreaseAllowed(double newf, double oldf, int nOuterIters) const;

  //! sets the output parameters for the linesearch in an "Output" sublist of the "Line Search" sublist. 
  bool setOutputParameters();  
  
protected:

  //! Algorithms used to determine convergence of the line search
  enum ConvergenceCriteriaType {
    //! Sufficient decrease condition
    ArmijoGoldstein, 
    //! Ared/Pred condition
    AredPred,
    //! Just accept the first step
    None
  };

  //! Interpolation Type
  enum InterpolationType {
    //! Use quadratic interpolation throughout
    Quadratic, 
    //! Use quadratic interpolation in the first inner iteration and cubic interpolation otherwise
    Cubic
  };

  //! Choice of convergence criteria
  ConvergenceCriteriaType convCriteria;

  //! Choice of interpolation type
  InterpolationType interpolationType;

  //! True is we should do at least one interpolation step
  bool doForceInterpolation;

  //! True is we should allow residual increases constrained by relIncrease
  bool allowIncrease;

  //! Allowable relative increase for residual norm
  double relIncrease;

  //! Maximum number relative increase steps allowed
  int numAllowed;

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
  NOX::Parameter::List* paramsPtr;

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
