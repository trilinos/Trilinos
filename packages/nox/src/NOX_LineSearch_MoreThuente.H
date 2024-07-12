// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_MORETHUENTE_H
#define NOX_LINESEARCH_MORETHUENTE_H

#include "NOX_LineSearch_Generic.H" // base class

#include "NOX_LineSearch_Utils_Printing.H"  // data class member
#include "NOX_LineSearch_Utils_Slope.H"     // data class member
#include "Teuchos_RCP.hpp"                  // data class member

// Forward Declarations
namespace NOX {
  namespace MeritFunction {
    class Generic;
  }
  class LineSearchCounters;
}

namespace NOX {
namespace LineSearch {

//! More'-Thuente Line Search. Original code by Dianne O'Leary, modfified by Tammy Kolda and Roger Pawlowski for the NOX project.  This version has been slightly optimized and also supports Homer Walker's work on adaptive forcing terms and Ared/Pred conditions.  It also allows for arbitrary merit functions and norms to be supplied by the user.
/*!

  This code is based on the More'-Thuente line search from the 1983
  MINPACK Project. More specifically, this code is based on Dianne
  O'Leary's 1991 Matlab-implementation of the More'-Thuente line
  search. The original comments are preserved in the descriptions of
  the individual subroutines. What follows is an updated summary.

  The merit function we are minimizing is given by

  \f[
  f(x) = 0.5 \|F(x)\|^2
  \f]

  (alternatively the user can define this)

  The purpose of the More'-Thuente line search is to find a step which
  satisfies a sufficient decrease condition and a curvature condition.
  At each stage the subroutine updates an interval of uncertainty with
  endpoints stx and sty. The interval of uncertainty is initially
  chosen so that it contains a minimizer of the modified function

  \f[
  f(x+{\rm stp} \; s) - f(x) - {\rm ftol} \; {\rm stp} \; (\nabla f(x)^T s).
  \f]

  If a step is obtained for which the modified function has a
  nonpositive function value and nonnegative derivative, then the
  interval of uncertainty is chosen so that it contains a minimizer of
  \f$f(x+{\rm stp}\;s)\f$.

  The algorithm is designed to find a step which satisfies one of two
  sufficient decrease conditions:

  (1) Armijo-Goldstein Condition
  \f[
  f(x + {\rm stp} \; s) \le f(x) + {\rm ftol} \; {\rm stp} \; (\nabla f(x)^T s),
  \f]

  or

  (2) Ared/Pred Condtition
  \f[
  F(x_{n-1}+ \lambda s) \le  F(x_{n-1})(1-\alpha(1-\eta))
  \f]


  and the curvature condition

  \f[
  \vert \nabla f(x + {\rm stp} \; s)^T s) \vert \le {\rm gtol} \;  |\nabla f(x)^T s|
  \f]

  If ftol is less than gtol and if, for example, the function is
  bounded below, then there is always a step which satisfies both
  conditions. If no step can be found which satisfies both conditions,
  then the algorithm usually stops when rounding errors prevent
  further progress. In this case stp only satisfies the sufficient
  decrease condition.

  <B> Modifications from NOX::LineSearch::MoreThuente </B> <br>
  1. Added the option to use Ared/Pred conditions as describe in Homer Walker's papers.
  2. Added support to use an adjustable forcing term as describe in Homer Walker's papers.
  3. Added the option to use directional derivatives in computing the slope instead of explicitly computing the Jacobian.  This eliminates the need to recompute the Jacobian at each inner iteration of the More'-Thuente algorithm.
  4. Added the ability to use the NOX::Parameter::UserNorm and NOX::Parameter::MeritFunction objects to supply user defined norms and merit functions to the line search.

  <B> Implementation </B> <br>
  This line search can be called via NOX::LineSearch::Manager.

  This line search is used if "More'-Thuente" is the "Method" in the
  "Line Search" sublist. (See NOX::LineSearch::Manager for details.)

  The following parameters can be specified for this line search in
  the "More'-Thuente" sublist of the "Line Search" sublist:

  - "Sufficient Decrease Condition" - Choice to use for the sufficient decrease condition. Options are "Ared/Pred" or "Armijo-Goldstein" (defaults to "Armijo-Goldstein").<br>
  1. "Armijo-Goldstein" conditions: \f$ f(x_{n-1}+ \lambda s) \le f(x_{n-1}) +\alpha \lambda f'(x_{n-1}) \f$ <br>
  2. "Ared/Pred" conditions: \f$ \| F(x_{n-1}+ \lambda s) \| \le \| F(x_{n-1}) \| (1-\alpha(1-\eta)) \f$ where \f$ \eta \f$ is the linear solve tolerance in the inexact Newton method. <br><br>

  - "Sufficient Decrease" - The ftol in the sufficient decrease
    condition (defaults to 1.0e-4)

  - "Curvature Condition" - The gtol in the curvature condition
    (defaults to 0.9999)

  - "Optimize Slope Calculation" - Boolean value.  If set to true the value of \f$ s^TJ^TF \f$ is estimated using a directional derivative in a call to NOX::LineSearch::Common::computeSlopeWithOutJac.  If false the slope computation is computed with the NOX::LineSearch::Common::computeSlope method.  Setting this to true eliminates having to compute the Jacobian at each inner iteration of the More'-Thuente line search (defaults to false).

  - "User Defined Norm" - The user can redefine the norm that is used in the
    Ared/Pred sufficient decrease condition by supplying a
    NOX::Parameter::UserNorm derived object in the parameter list with
    this key.

  - "Merit Function" - The user can supply their own merit function to the
    line search by supplying a NOX::Parameter::MeritFunction derived object
    with this key.

  - "Interval Width" - The maximum width of the interval containing
    the minimum of the modified function (defaults to 1.0e-15)

  - "Maximum Step" - maximum allowable step length (defaults to 1.0e6)

  - "Minimum Step" - minimum allowable step length (defaults to 1.0e-12)

  - "Max Iters" - maximum number of right-hand-side and corresponding
    Jacobian evaluations (defaults to 20)

  - "Default Step" - starting step length (defaults to 1.0)

  - "Recovery Step Type" - Determines the step size to take when the
    line search fails.  Choices are:

     <ul>
     <li> "Constant" [default] - Uses a constant value set in "Recovery Step".
     <li> "Last Computed Step" - Uses the last value computed by the
                                 line search algorithm.
     </ul>

  - "Recovery Step" - The value of the step to take when the line
    search fails. Only used if the "Recovery Step Type" is set to
    "Constant". Defaults to value for "Default Step".

  <B> Output Parameters </B> <br>
A sublist for output parameters will be created called "Output" in the parameter list used to instantiate or reset the class.  Valid output parameters are:

  - "Total Number of Line Search Calls" - Total number of calls to the compute() method of this line search.

  - "Total Number of Non-trivial Line Searches" - The total number of steps that could not directly take a full step and meet the required "Convergence Criteria" (i.e. The line search had to reduce the step length using inner iteration calculations over iterate \f$ k \f$).

  - "Total Number of Failed Line Searches" - total number of line searches that failed and used a recovery step.

  - "Total Number of Line Search Inner Iterations" - total number of inner iterations \f$ k \f$ performed by this object.



*/

class MoreThuente : public Generic {

public:

  //! Constructor
  MoreThuente(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params);

  //! Destructor
  ~MoreThuente();

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
           const NOX::Abstract::Vector& dir,
           const NOX::Solver::Generic& s);


private:

  //! More'-Thuente's cvsrch subroutine.
  /*!
    This translation into C++ is based on a translation into Matlab by Dianne O'Leary.
    The return code corresponds to the "info" code described below.
    The original documentation follows.

\verbatim
%   Translation of minpack subroutine cvsrch
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cvsrch
%
%     The purpose of cvsrch is to find a step which satisfies
%     a sufficient decrease condition and a curvature condition.
%     The user must provide a subroutine which calculates the
%     function and the gradient.
%
%     At each stage the subroutine updates an interval of
%     uncertainty with endpoints stx and sty. The interval of
%     uncertainty is initially chosen so that it contains a
%     minimizer of the modified function
%
%          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
%
%     If a step is obtained for which the modified function
%     has a nonpositive function value and nonnegative derivative,
%     then the interval of uncertainty is chosen so that it
%     contains a minimizer of f(x+stp*s).
%
%     The algorithm is designed to find a step which satisfies
%     the sufficient decrease condition
%
%           f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),
%
%     and the curvature condition
%
%           abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).
%
%     If ftol is less than gtol and if, for example, the function
%     is bounded below, then there is always a step which satisfies
%     both conditions. If no step can be found which satisfies both
%     conditions, then the algorithm usually stops when rounding
%     errors prevent further progress. In this case stp only
%     satisfies the sufficient decrease condition.
%
%     The subroutine statement is
%
%        subroutine cvsrch(fcn,n,x,f,g,s,stp,ftol,gtol,xtol,
%                          stpmin,stpmax,maxfev,info,nfev,wa)
%     where
%
%    fcn is the name of the user-supplied subroutine which
%         calculates the function and the gradient.  fcn must
%            be declared in an external statement in the user
%         calling program, and should be written as follows.
%
%      subroutine fcn(n,x,f,g)
%         integer n
%         f
%         x(n),g(n)
%      ----------
%         Calculate the function at x and
%         return this value in the variable f.
%         Calculate the gradient at x and
%         return this vector in g.
%      ----------
%      return
%      end
%
%       n is a positive integer input variable set to the number
%      of variables.
%
%    x is an array of length n. On input it must contain the
%      base point for the line search. On output it contains
%         x + stp*s.
%
%    f is a variable. On input it must contain the value of f
%         at x. On output it contains the value of f at x + stp*s.
%
%    g is an array of length n. On input it must contain the
%         gradient of f at x. On output it contains the gradient
%         of f at x + stp*s.
%
%    s is an input array of length n which specifies the
%         search direction.
%
%    stp is a nonnegative variable. On input stp contains an
%         initial estimate of a satisfactory step. On output
%         stp contains the final estimate.
%
%       ftol and gtol are nonnegative input variables. Termination
%         occurs when the sufficient decrease condition and the
%         directional derivative condition are satisfied.
%
%    xtol is a nonnegative input variable. Termination occurs
%         when the relative width of the interval of uncertainty
%      is at most xtol.
%
%    stpmin and stpmax are nonnegative input variables which
%      specify lower and upper bounds for the step.
%
%    maxfev is a positive integer input variable. Termination
%         occurs when the number of calls to fcn is at least
%         maxfev by the end of an iteration.
%
%    info is an integer output variable set as follows:
%
%      info = 0  Improper input parameters.
%
%      info = 1  The sufficient decrease condition and the
%                   directional derivative condition hold.
%
%      info = 2  Relative width of the interval of uncertainty
%            is at most xtol.
%
%      info = 3  Number of calls to fcn has reached maxfev.
%
%      info = 4  The step is at the lower bound stpmin.
%
%      info = 5  The step is at the upper bound stpmax.
%
%      info = 6  Rounding errors prevent further progress.
%                   There may not be a step which satisfies the
%                   sufficient decrease and curvature conditions.
%                   Tolerances may be too small.
%
%       nfev is an integer output variable set to the number of
%         calls to fcn.
%
%    wa is a work array of length n.
%
%     Subprograms called
%
%    user-supplied......fcn
%
%    MINPACK-supplied...cstep
%
%    FORTRAN-supplied...abs,max,min
%
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
\endverbatim


   */
  int cvsrch(NOX::Abstract::Group& newgrp, double& step,
         const NOX::Abstract::Group& oldgrp,
         const NOX::Abstract::Vector& dir,
         const Solver::Generic& s);


  //! More'-Thuente's cstep subroutine.
  /*!
    \verbatim
%   Translation of minpack subroutine cstep
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cstep
%
%     The purpose of cstep is to compute a safeguarded step for
%     a linesearch and to update an interval of uncertainty for
%     a minimizer of the function.
%
%     The parameter stx contains the step with the least function
%     value. The parameter stp contains the current step. It is
%     assumed that the derivative at stx is negative in the
%     direction of the step. If brackt is set true then a
%     minimizer has been bracketed in an interval of uncertainty
%     with endpoints stx and sty.
%
%     The subroutine statement is
%
%       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
%                        stpmin,stpmax,info)
%
%     where
%
%       stx, fx, and dx are variables which specify the step,
%         the function, and the derivative at the best step obtained
%         so far. The derivative must be negative in the direction
%         of the step, that is, dx and stp-stx must have opposite
%         signs. On output these parameters are updated appropriately.
%
%       sty, fy, and dy are variables which specify the step,
%         the function, and the derivative at the other endpoint of
%         the interval of uncertainty. On output these parameters are
%         updated appropriately.
%
%       stp, fp, and dp are variables which specify the step,
%         the function, and the derivative at the current step.
%         If brackt is set true then on input stp must be
%         between stx and sty. On output stp is set to the new step.
%
%       brackt is a logical variable which specifies if a minimizer
%         has been bracketed. If the minimizer has not been bracketed
%         then on input brackt must be set false. If the minimizer
%         is bracketed then on output brackt is set true.
%
%       stpmin and stpmax are input variables which specify lower
%         and upper bounds for the step.
%
%       info is an integer output variable set as follows:
%         If info = 1,2,3,4,5, then the step has been computed
%         according to one of the five cases below. Otherwise
%         info = 0, and this indicates improper input parameters.
%
%     Subprograms called
%
%       FORTRAN-supplied ... abs,max,min,sqrt
%                        ... dble
%
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
    \endverbatim
  */
  int cstep(double& stx, double& fx, double& dx,
        double& sty, double& fy, double& dy,
        double& stp, double& fp, double& dp,
        bool& brackt, double stmin, double stmax);

  //! max
  double max(double a, double b);

  //! min
  double min(double a, double b);

  //! absmax - returns the max of the absolute value of the three inputs
  double absmax(double a, double b, double c);


private:

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Common line search printing utilities.
  NOX::LineSearch::Utils::Printing print;

  //! Common common counters for line searches.
  NOX::LineSearchCounters* counter;

  //! Common slope calculations for line searches.
  NOX::LineSearch::Utils::Slope slope;

  //! Pointer to the input parameter list.  This is saved so data stored in the counter object can be written to the "Output" parameter list.
  Teuchos::ParameterList* paramsPtr;

  //! Sufficient decrease parameter
  double ftol;
  //! Curvature condition parameter
  double gtol;
  //! Threshold on width of the uncertainty
  double xtol;
  //! Minimum step
  double stpmin;
  //! Maximum step
  double stpmax;
  //! Maximum number of function evaluations
  int maxfev;

  //! Type of recovery step to use
  enum RecoveryStepType {
    //! Use a constant value
    Constant,
    //! Use the last value computed in the line search algorithm
    LastComputedStep
  };

  //! Choice of the recovery step type; uses "Recovery Step Type" parameter
  RecoveryStepType recoveryStepType;

  //! Default step length
  double defaultstep;
  //! Recovery step if linesearch fails
  double recoverystep;

  //! Options for sufficient decrease criteria
  enum SufficientDecreaseCondition {ArmijoGoldstein, AredPred};

  //! Chioce for the sufficient decrease criteria
  SufficientDecreaseCondition suffDecrCond;

  //! If set to true the value of \f$ s^TJ^Tf \f$ is estimated using a directional derivative.  This eliminates having to compute the Jacobian at each inner iteration of the line search.
  bool useOptimizedSlopeCalc;

  //! Stores the merit function.
  Teuchos::RCP<NOX::MeritFunction::Generic> meritFuncPtr;

};
} // namespace LineSearch
} // namespace NOX
#endif
