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

#ifdef WITH_PRERELEASE

#ifndef NOX_LINESEARCH_POLYNOMIALNEW_H
#define NOX_LINESEARCH_POLYNOMIALNEW_H

#include "NOX_LineSearch_Generic.H" // base class

#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member

namespace NOX {
namespace Parameter {
class UserNorm;
class MeritFunction;
}
}

namespace NOX {
namespace LineSearch {

/*!
  \brief
  A polynomial line search, either quadratic or cubic.

  This line search can be called via NOX::LineSearch::Manager.
 
  The goal of the line search is to find a step \f$\lambda\f$ for the
  calculation \f$x_{\rm new} = x_{\rm old} + \lambda d\f$, given
  \f$x_{\rm old}\f$ and \f$ d \f$.  To accomplish this goal, we
  minimize a merit function \f$ \phi(\lambda) \f$ that measures the
  "goodness" of the step \f$\lambda\f$.  The standard merit function is

  \f[ 
  \phi(\lambda) \equiv \frac{1}{2}||F (x_{\rm old} + \lambda s)||^2,
  \f]

  but a user defined merit function can be used instead (see
  NOX::Parameter::MeritFunction for details).  Our first attempt is
  always to try a step \f$ \lambda_0 \f$, and then check the stopping
  criteria. The value of \f$ \lambda_0 \f$ is the specified by the
  "Default Step" parameter. If the first try doesn't work, then we
  successively minimize polynomial approximations, \f$ p_k(\lambda)
  \approx \phi(\lambda) \f$.

  <b>Stopping Criteria</b>

  The inner iterations continue until:

  <ul>

  <li>The sufficient decrease condition is met; see checkConvergence()
  for details.

  <li>The maximum iterations are reached; see parameter 
  "Max Iters". 
  This is considered a failure and the recovery step is
  used; see parameter "Recovery Step".

  <li>The minimum step length is reached; see parameter 
  "Minimum Step". 
  This is considered a line search failure
  and the recovery step is used; see parameter 
  "Recovery Step".

  </ul>

  <b> %Polynomial Models of the Merit Function </b>

  We compute \f$ p_k(\lambda) \f$ by interpolating \f$ \phi(0) \f$ at
  \f$ \phi'(0) \f$ and \f$ \phi(\lambda_{k-1}) \f$ where \f$ \lambda_{k-1} \f$
  is the \f$ k-1 \f$st approximation to the step. 

  The first interpolation step is always quadratic.  Depending on the
  choice of "Interpolation Type" parameter, the
  following steps are either quadratic or cubic.
  
  If we are using a <b>quadratic model</b>, the \f$k\f$-th step approximation is
  calculated as

  \f[ 
  \lambda_{k} = \frac{-\phi'(0) \lambda_{k-1}^2 }{2 \left[ \phi(\lambda_{k-1}) - \phi(0) 
  -\phi'(0) \lambda_{k-1} \right]} 
  \f]

  See computeSlope() for details on how \f$ \phi'(0) \f$ is
  calculated.

  If we are using a <b>cubic model</b>, the \f$k\f$-th step (\f$ k>1
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

  <B> Bounds on the step length </B> 
 
  We do not allow the step to grow or shrink too quickly, using the following bounds:

  \f[ 
  \gamma_{min} \; \lambda_{k-1} \leq \lambda_k \le \gamma_{max} \; \lambda_{k-1}
  \f]

  Here \f$ \gamma_{min} \f$ and \f$ \gamma_{max} \f$ are defined by
  parameters "Min Bounds Factor" and 
  "Max Bounds Factor".

  <B> Input Parameters </B> 

  "Line Search"/"PolynomialNew":

  <ul>
 
  <li> \anchor defaultstep "Default Step" - starting step length (defaults to 1.0)

  <li> \anchor minimumstep "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)

  <li> \anchor recovertstep "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")

  <li> \anchor maxiters "Max Iters" - maximum number of line search iterations (default to 100)

  <li> \anchor sufficientdecreasecondition "Sufficient Decrease
    Condition" - choice of algorithm used to determine if the line
    search algorithm has converged. Choices are (See above for
    explanations):
    

    <ol>

    <li> "Armijo-Goldstein" [default]
    <li> "Ared/Pred" 
    <li> "None"

    </ol>

  <li> "User Defined Norm" - The user can redefine the norm that is
    used in the Ared/Pred sufficient decrease condition; see
    NOX::Parameter::UserNorm for details.

  <li> "User Defined Merit Function" - The user can redefine the merit
  function used; see NOX::Parameter::MeritFunction for details.

  <li> \anchor alphafactor "Alpha Factor" - Used in Armijo-Goldstein stopping criteria \f$
    \alpha \f$ (defaults to 1.0e-4)

  <li> \anchor minboundsfactor "Min Bounds Factor" - \f$ \gamma_{min} \f$, the factor that limits
    the minimum size of the new step based on the previous step
    (defaults to 0.1)

  <li> \anchor maxboundsfactor "Max Bounds Factor" - \f$ \gamma_{max} \f$, the factor that limits
    the maximum size of the new step based on the previous step
    (defaults to 0.5)

  <li> \anchor forceinterpolation "Force Interpolation" (boolean) -
    Set to true at least one interpolation step should be used. The
    default is false which means that the line search will stop if the
    default step length satisfies the convergence criteria. (defaults
    to false)

  <li> \anchor interpolationtype "Interpolation Type" - Type of
    interpolation that should be used. Choices are
    
     <ol>
     <li> "Cubic" 
     <li> "Quadratic"
     </ol>

  <li> \anchor usecounters "Use Counters" (boolean) - True if we
  should use counters and then output the result to the paramter list
  as described in \ref outputparameters "Output Parameters". Defaults
  to true.

  <li> "Allowed Relative Increase" - See #maxRelativeIncrease and
  checkConvergence().  Defaults to 100.

  <li> "Maximum Iteration for Increase" - Maximum index of the
  nonlinear iteration for which we allow a relative increase. See
  #maxIncreaseIter and checkConvergence(). Defaults to zero.

  </ul>   

  \anchor outputparameters <B> Output Parameters </B> 

  If the "Use Counters" parameter is set to true, then a sublist
  for output parameters called "Output" will be created in the
  parameter list used to instantiate or reset the class.  Valid output
  parameters are:

  <ul>

  <li> "Total Number of Line Search Calls" - Total number of calls to
  the compute() method of this line search.

  <li> "Total Number of Non-trivial Line Searches" - The total number of
    steps that could not directly take a full step and meet the
    required "Sufficient Decrease Condition" (i.e. The line search had to
    reduce the step length using inner iteration calculations over
    iterate \f$ k \f$).

  <li> "Total Number of Failed Line Searches" - total number of line
    searches that failed and used a recovery step.

  <li> "Total Number of Line Search Inner Iterations" - total number of
    inner iterations \f$ k \f$ performed by this object.

  </ul>

  <b>References</b>

  This line search is based on materials in the following:

  - Section 8.3.1 in C.T. Kelley, "Iterative Methods for Linear and
    Nonlinear Equations", SIAM, 1995.

  - Section 6.3.2 and Algorithm 6.3.1 of J. E. Dennis Jr. and Robert
    B.  Schnabel, "Numerical Methods for Unconstrained Optimization
    and Nonlinear Equations," Prentice Hall, 1983.

  - Section 3.4 of Jorge Nocedal and Stephen J. Wright, "Numerical
    Optimization,"Springer, 1999.

  - "An Inexact Newton Method for Fully Coupled Solution of the Navier-Stokes 
    Equations with Heat and Mass Transfer", Shadid, J. N., Tuminaro, R. S., 
    and Walker, H. F., Journal of Computational Physics, 137, 155-185 (1997)

  \author Russ Hooper, Roger Pawlowski, Tammy Kolda
*/

class PolynomialNew : public Generic {

public:

  //! Constructor
  PolynomialNew(const NOX::Utils& u, NOX::Parameter::List& params);

  //! Destructor
  ~PolynomialNew();
  
  // derived
  bool reset(NOX::Parameter::List& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step, 
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s);
	       
protected:

  //! Returns true if converged.
  /*!  

    If the "Force Interpolation" parameter is true and this is the
    first nonlinear iteration (i.e., nNonlinearIterations is 1), then
    returns false.

    Next, it checks to see if the "Relative Increase" condition is
    satisfied. If so, then it returns true (converged).

    The "Relative Increase" condition is satisfied if the following
    holds.
    -# #doAllowIncrease is true
    -# nNonlinearIters <= #maxRelativeIncrease
    -# (newValue / oldValue) < #maxRelativeIncrease

    Finally, it checks the actual convergence conditions.
    
    In the case of <b>Armijo-Goldstein</b> criteria, we return true if
    \f[
    \phi(\lambda) \leq \phi(0) + \alpha \; \lambda \; \phi'(0)
    \f]
    where 
    - newValue = \f$ \phi(\lambda)  \f$
    - oldValue = \f$ \phi(0) \f$
    - oldSlope = \f$ \phi'(0) \f$
    - step = \f$ \lambda  \f$

    In the case of <b>Ared/Pred</b>, we return true if
    \f[
    \| F(x_{\rm old} + \lambda d )\| \leq \| F(x_{\rm old}) \| (1 - \alpha (1 - \eta))
    \f]
    where
    - newValue = \f$ \| F(x_{\rm old} + \lambda d )\|  \f$
    - oldValue = \f$ \| F(x_{\rm old}) \| \f$
    - eta = Value of \f$\eta\f$ for last forcing term calculation in NOX::Direction::Newton
    - step = \f$ \lambda \f$
  */
  bool checkConvergence(double newValue, double oldValue, double oldSlope, 
			double step, double eta, int nNonlinearIters) const;

  //! Updates the newGrpand computes the values of newf_1, newf_2 and newf_interp
  bool updateGrp(NOX::Abstract::Group& newGrp, 
		 const NOX::Abstract::Group& oldGrp,
		 const NOX::Abstract::Vector& dir, 
		 double step) const;


  //! Computes the value of the merit function
  /*!
    If a "User Defined Merit Function" parameter is specified, then
    the NOX::Parameter::MeritFunction::computeF function on the
    user-supplied merit function is used.

    Otherwise, we compute and return

    \f[
    \frac{1}{2} \| F(x) \|^2
    \f]
   */
  double computePhi(const NOX::Abstract::Group& grp);

  //! Compute the value used to determine sufficient decrease
  /*!
    If the "Sufficient Decrease Condition" is set to "Ared/Pred", then
    we do the following.  If a "User Defined Norm" parameter is
    specified, then the NOX::Parameter::UserNorm::norm function on the
    user-supplied merit function is used. If the user does not provide
    a norm, we return \f$ \|F(x)\|  \f$.

    If the "Sufficient Decrease Condition" is <em>not</em> set to
    "Ared/Pred", then we simply return phi.

    \param phi - Should be equal to computePhi(grp).
  */
  double computeValue(const NOX::Abstract::Group& grp, double phi);


  //! Compute \f$ \phi'(x) \f$
  /*!

    If a "User Defined Merit Function" is provided, then we use
    NOX::Parameter::MeritFunction::computeSlope to compute this value.

    Otherwise, we are using the standard merit function, so we want to
    compute
    \f[
    d^T (J(x)^T F(x))
    \f]
    To do this, we use NOX::LineSearch::Utils::Slope.

   */
  double computeSlope(const NOX::Abstract::Vector& dir, 
		      const NOX::Abstract::Group& grp);

  //! Used to print opening remarks for each call to compute().
  void printOpeningRemarks() const;

  //! Prints a warning message saying that the slope is negative
  void printBadSlopeWarning(double slope) const;

protected:

  //! Algorithms used to determine termination of the line search
  enum SufficientDecreaseType {
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

  //! Choice for sufficient decrease condition
  SufficientDecreaseType suffDecrCond;

  //! Choice of interpolation type
  InterpolationType interpolationType;

   //! Minimum step length (i.e., when we give up)
  double minStep;
 
 //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  //! Maximum iterations
  int maxIters;

  //! The \f$ \alpha \f$ for the Armijo-Goldstein condition, or the
  //! \f$ \alpha \f$ for the Ared/Pred condition.
  double alpha;

  //! Factor that limits the minimum size of the new step based on the previous step.
  double minBoundFactor;

  //! Factor that limits the maximum size of the new step based on the previous step.
  double maxBoundFactor;

  //! True is we should force at least one interpolation step
  bool doForceInterpolation;

  //! True if we sometimes allow an increase(!) in the decrease measure
  bool doAllowIncrease;

  //! Maximum allowable relative increase for decrease meausre (if allowIncrease is true)
  double maxRelativeIncrease;

  //! No increases are allowed if the number of nonlinear iterations is greater than this value
  int maxIncreaseIter;

  //! True if we should use #counterUtil and output the results
  bool useCounter;

  //! Pointer to the input parameter list. 
  /*! We need this to later create an "Output" sublist to store output
      parameters from #counterUtil. 

      This should not be deleted in the destructor 
  */
  NOX::Parameter::List* paramsPtr;

  //! Common line search printing utilities.
  NOX::LineSearch::Utils::Printing print;

  //! Common common counters for line searches.
  NOX::LineSearch::Utils::Counters counter;

  //! Common slope calculations for line searches.
  NOX::LineSearch::Utils::Slope slopeUtil;

  //! Pointer to a user supplied norm.
  /*! This should not be deleted in the destructor */
  NOX::Parameter::UserNorm* userNormPtr;

  //! Pointer to a user supplied merit function.
  /*! This should not be deleted in the destructor */
  NOX::Parameter::MeritFunction* meritFuncPtr;

};
} // namespace LineSearch
} // namespace NOX
#endif
#endif
