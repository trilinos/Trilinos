// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_POLYNOMIALNEW_H
#define NOX_LINESEARCH_POLYNOMIALNEW_H

#include "NOX_LineSearch_Generic.H" // base class

#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member
#include "Teuchos_RCP.hpp"                  // class data member

// Forward Declarations
namespace NOX {
  namespace MeritFunction {
    class Generic;
  }
  class LineSearchCounters;
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
  computePhi() for details).  Our first attempt is
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

  We compute \f$ p_k(\lambda) \f$ by interpolating \f$f\f$. For the
  quadratic fit, we interpolate \f$ \phi(0) \f$, \f$ \phi'(0) \f$, and
  \f$ \phi(\lambda_{k-1}) \f$ where \f$ \lambda_{k-1} \f$ is the \f$
  k-1 \f$st approximation to the step. For the cubit fit, we
  additionally include \f$\phi(\lambda{k-2})\f$.

  The steps are calculated iteratively as follows, depending on the
  choice of "Interpolation Type".

  <ul>

  <li> "Quadratic" - We construct a quadratic model of \f$\phi\f$, and solve for \f$\lambda\f$ to get

  \f[
  \lambda_{k} = \frac{-\phi'(0) \lambda_{k-1}^2 }{2 \left[ \phi(\lambda_{k-1}) - \phi(0)
  -\phi'(0) \lambda_{k-1} \right]}
  \f]

  <li> "Cubic" - We construct a cubic model of \f$\phi\f$, and solve for
  \f$\lambda\f$. We use the quadratic model to solve for \f$\lambda_1\f$; otherwise,

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

  <li> "Quadratic3" - We construct a quadratic model of \f$\phi\f$ using
  \f$\phi(0)\f$, \f$ \phi(\lambda_{k-1}) \f$ , and
  \f$\phi(\lambda_{k-2})\f$. No derivative information for \f$\phi\f$
  is used. We let \f$\lambda_1 = \frac{1}{2} \lambda_0\f$, and otherwise

  \f[
  \lambda_k = - \frac{1}{2}
  \frac{\lambda_{k-1}^2 [\phi(\lambda_{k-2}) -\phi(0)]
  - \lambda_{k-2}^2 [\phi(\lambda_{k-1}) -\phi(0)]}
  {\lambda_{k-2} [\phi(\lambda_{k-1}) -\phi(0)]
  - \lambda_{k-1} [\phi(\lambda_{k-2}) -\phi(0)]}
  \f]

  </ul>

  For "Quadratic" and "Cubic", see computeSlope() for details on how
  \f$ \phi'(0) \f$ is calculated.

  <B> Bounds on the step length </B>

  We do not allow the step to grow or shrink too quickly by enforcing
  the following bounds:

  \f[
  \gamma_{min} \; \lambda_{k-1} \leq \lambda_k \le \gamma_{max} \; \lambda_{k-1}
  \f]

  Here \f$ \gamma_{min} \f$ and \f$ \gamma_{max} \f$ are defined by
  parameters "Min Bounds Factor" and
  "Max Bounds Factor".

  <B> Input Parameters </B>

  "Line Search":

  <ul>
  <li>"Method" = "Polynomial" [required]
  </ul>

  "Line Search"/"Polynomial":

  <ul>

  <li> "Default Step" - Starting step length, i.e., \f$\lambda_0\f$.
  Defaults to 1.0.

  <li> "Max Iters" - Maximum number of line search iterations. The
  search fails if the number of iterations exceeds this
  value. Defaults to 100.

  <li> "Minimum Step" - Minimum acceptable step length. The search
  fails if the computed \f$\lambda_k\f$ is less than this
  value. Defaults to 1.0e-12.

  <li> "Recovery Step Type" - Determines the step size to take when the
    line search fails.  Choices are:

     <ul>
     <li> "Constant" [default] - Uses a constant value set in "Recovery Step".
     <li> "Last Computed Step" - Uses the last value computed by the
                                 line search algorithm.
     </ul>

  <li> "Recovery Step" - The value of the step to take when the line
    search fails. Only used if the "Recovery Step Type" is set to
    "Constant". Defaults to value for "Default Step".

  <li> "Interpolation Type" - Type of interpolation that should be
    used. Choices are

     <ul>
     <li> "Cubic" [default]
     <li> "Quadratic"
     <li> "Quadratic3"
     </ul>

  <li> "Min Bounds Factor" - Choice for \f$ \gamma_{min} \f$, i.e.,
    the factor that limits the minimum size of the new step based on
    the previous step. Defaults to 0.1.

  <li> "Max Bounds Factor" - Choice for \f$ \gamma_{max} \f$, i.e.,
    the factor that limits the maximum size of the new step based on
    the previous step.  Defaults to 0.5.

  <li> "Sufficient Decrease Condition" - See checkConvergence() for
  details. Choices are:

    <ul>

    <li> "Armijo-Goldstein" [default]
    <li> "Ared/Pred"
    <li> "None"

    </ul>

  <li> "Alpha Factor" - %Parameter choice for sufficient decrease
  condition. See checkConvergence() for details. Defaults to 1.0e-4.

  <li> "Force Interpolation" (boolean) - Set to true if at least one
    interpolation step should be used. The default is false which
    means that the line search will stop if the default step length
    satisfies the convergence criteria. Defaults to false.

  <li> "Use Counters" (boolean) - Set to true if we should use
  counters and then output the result to the paramter list as
  described in \ref outputparameters "Output Parameters". Defaults to
  true.

  <li> "Maximum Iteration for Increase" - Maximum index of the
  nonlinear iteration for which we allow a relative increase. See
  checkConvergence() for further details. Defaults to 0 (zero).

  <li> "Allowed Relative Increase" - See checkConvergence() for
  details.  Defaults to 100.

  <li> "User Defined Merit Function" - The user can redefine the merit
  function used; see computePhi() and NOX::Parameter::MeritFunction
  for details.

  <li> "User Defined Norm" - The user can redefine the norm that is
    used in the Ared/Pred sufficient decrease condition; see
    computeValue() and NOX::Parameter::UserNorm for details.

  </ul>

  \anchor outputparameters <B> Output Parameters </B>

  If the "Use Counters" parameter is set to true, then a sublist
  for output parameters called "Output" will be created in the
  parameter list used to instantiate or reset the class.  Valid output
  parameters are:

  <ul>

  <li> "Total Number of Line Search Calls" - Total number of calls to
  the compute() method of this line search.

  <li> "Total Number of Non-trivial Line Searches" - Total number of
    steps that could not directly take a full step and meet the
    required "Sufficient Decrease Condition" (i.e., the line search had to
    do at least one interpolation).

  <li> "Total Number of Failed Line Searches" - Total number of line
    searches that failed and used a recovery step.

  <li> "Total Number of Line Search Inner Iterations" - Total number of
    inner iterations of all calls to compute().

  </ul>

  <b>References</b>

  This line search is based on materials in the following:

  - Section 8.3.1 in C.T. Kelley, "Iterative Methods for %Linear and
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

class Polynomial : public Generic {

public:

  //! Constructor
  Polynomial(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  //! Destructor
  ~Polynomial();

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
           const NOX::Abstract::Vector& dir,
           const NOX::Solver::Generic& s);

protected:

  //! Returns true if converged.
  /*!
    We go through the following checks for convergence.

    <ol>

    <li> If the "Force Interpolation" parameter is true and this is the
    first inner iteration (i.e., nIters is 1), then we return \b false.

    <li> Next, it checks to see if the "Relative Increase" condition
    is satisfied. If so, then we return \b true.  The "Relative
    Increase" condition is satisfied if \e both of the following two
    conditions hold.

    - The number of nonlinear iterations is less than or equal to the
       value specified in the "Maximum Iteration for Increase"
       parameter

    - The ratio of newValue to oldValue is less than the value
       specified in "Allowed Relative Increase".

    <li> Last, we check the sufficient decrease condition. We return
    \b true if it's satisfied, and \b false otherwise. The condition
    depends on the value of "Sufficient Decrease Condition",
    as follows.

    <ul>
    <li> "Armijo-Goldstein": Return true if
    \f[
    \phi(\lambda) \leq \phi(0) + \alpha \; \lambda \; \phi'(0)
    \f]

    <li>"Ared/Pred": Return true if
    \f[
    \| F(x_{\rm old} + \lambda d )\| \leq \| F(x_{\rm old}) \| (1 - \alpha (1 - \eta))
    \f]


    <li> "None": Always return true.
    </ul>

    For the first two cases, \f$\alpha\f$ is specified by the
    parameter "Alpha Factor".

    </ol>

    \param newValue Depends on the "Sufficient Decrease Condition" parameter.
    <ul>
    <li>"Armijo-Goldstein" - \f$ \phi(\lambda) \f$
    <li>"Ared/Pred" - \f$ \| F(x_{\rm old} + \lambda d )\| \f$
    <li>"None" - N/A
    </ul>

    \param oldValue Depends on the "Sufficient Decrease Condition" parameter.
    <ul>
    <li>"Armijo-Goldstein" - \f$ \phi(0) \f$
    <li>"Ared/Pred" - \f$ \| F(x_{\rm old} )\| \f$
    <li>"None" - N/A
    </ul>

    \param oldSlope Only applicable to "Armijo-Goldstein", in which
    case it should be \f$\phi'(0)\f$.

    \param step The current step, \f$\lambda\f$.

    \param eta Only applicable to "Ared/Pred", in which case it
    should be the value of \f$\eta\f$ for last forcing term
    calculation in NOX::Direction::Newton

    \param nIters Number of inner iterations.

    \param nNonlinearIters Number of nonlinear iterations.

    \note The norm used for "Ared/Pred" can be specified by the user
    by using the "User Defined Norm" parameter; this parameter is any
    object derived from NOX::Parameter::UserNorm.

  */
  bool checkConvergence(double newValue, double oldValue, double oldSlope,
            double step, double eta, int nIters, int nNonlinearIters) const;

  //! Updates the newGrp by computing a new x and a new F(x)
  /*!
    Let

    - \f$x_{\rm new}\f$ denote the new solution to be calculated (corresponding to \c newGrp)
    - \f$x_{\rm old}\f$ denote the previous solution (i.e., the result of oldGrp.getX())
    - \f$d\f$ denotes the search direction (\c dir).
    - \f$\lambda\f$ denote the step (\c step),

    We compute
    \f[
    x_{\rm new} = x_{\rm old} + \lambda d,
    \f]
    and \f$ F(x_{\rm new})\f$. The results are stored in \c newGrp.
   */
  bool updateGrp(NOX::Abstract::Group& newGrp,
         const NOX::Abstract::Group& oldGrp,
         const NOX::Abstract::Vector& dir,
         double step) const;


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

  //! Used to print opening remarks for each call to compute().
  void printOpeningRemarks() const;

  //! Prints a warning message saying that the slope is negative
  void printBadSlopeWarning(double slope) const;

protected:

  //! Types of sufficient decrease conditions used by checkConvergence()
  enum SufficientDecreaseType {
    //! Armijo-Goldstein conditions
    ArmijoGoldstein,
    //! Ared/Pred condition
    AredPred,
    //! No condition
    None
  };

  //! Interpolation types used by compute().
  enum InterpolationType {

    //! Use quadratic interpolation throughout
    Quadratic,

    //! Use quadratic interpolation in the first inner iteration and cubic interpolation otherwise
    Cubic,

    //! Use a 3-point quadratic line search. (The second step is 0.5 times the default step.)
    Quadratic3
  };

  //! Type of recovery step to use
  enum RecoveryStepType {
    //! Use a constant value
    Constant,
    //! Use the last value computed in the line search algorithm
    LastComputedStep
  };

  //! Choice for sufficient decrease condition; uses "Sufficient Decrease Condition" parameter
  SufficientDecreaseType suffDecrCond;

  //! Choice of interpolation type; uses "Interpolation Type" parameter
  InterpolationType interpolationType;

  //! Choice of the recovery step type; uses "Recovery Step Type" parameter
  RecoveryStepType recoveryStepType;

   //! Minimum step length (i.e., when we give up); uses "Mimimum Step" parameter
  double minStep;

  //! Default step; uses "Default Step"  parameter
  double defaultStep;

  //! Default step for linesearch failure; uses "Recovery Step" parameter
  double recoveryStep;

  //! Maximum iterations; uses "Max Iters" parameter
  int maxIters;

  /*! \brief The \f$ \alpha \f$ for the Armijo-Goldstein condition, or
      the \f$ \alpha \f$ for the Ared/Pred condition; see
      checkConvergence().  Uses "Alpha Factor" parameter. */
  double alpha;

  /*! \brief Factor that limits the minimum size of the new step based
      on the previous step; uses "Min Bounds Factor" parameter */
  double minBoundFactor;

  /*! \brief Factor that limits the maximum size of the new step based
      on the previous step; uses "Max Bounds Factor" parameter. */
  double maxBoundFactor;

  /*! \brief True is we should force at least one interpolation step;
      uses "Force Interpolation" parameter. */
  bool doForceInterpolation;

  /*! \brief No increases are allowed if the number of nonlinear
    iterations is greater than this value; uses "Maximum Iterations
    for Increase" */
  int maxIncreaseIter;

  /*! \brief True if we sometimes allow an increase(!) in the decrease
    measure, i.e., maxIncreaseIter > 0. */
  bool doAllowIncrease;

  /*! \brief Maximum allowable relative increase for decrease meausre (if
    allowIncrease is true); uses "Allowed Relative Increase" parameter */
  double maxRelativeIncrease;

  /*! \brief True if we should use #counter and output the results; uses
    "Use Counters" parameter.*/
  bool useCounter;

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Pointer to the input parameter list.
  /*! We need this to later create an "Output" sublist to store output
      parameters from #counter.
  */
  Teuchos::ParameterList* paramsPtr;

  //! Common line search printing utilities.
  NOX::LineSearch::Utils::Printing print;

  //! Common common counters for line searches.
  NOX::LineSearchCounters* counter;

  //! Common slope calculations for line searches.
  NOX::LineSearch::Utils::Slope slopeUtil;

  //! Pointer to a user supplied merit function.
  Teuchos::RCP<NOX::MeritFunction::Generic> meritFuncPtr;

};
} // namespace LineSearch
} // namespace NOX
#endif
