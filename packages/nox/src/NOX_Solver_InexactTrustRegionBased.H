// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SOLVER_INEXACTTRUSTREGIONBASED_H
#define NOX_SOLVER_INEXACTTRUSTREGIONBASED_H

#include "NOX_Solver_Generic.H"        // base class

// Class data elements
#include "Teuchos_ParameterList.hpp"
#include "NOX_Direction_Utils_InexactNewton.H"
#include "NOX_SolverStats.hpp"
#include "Teuchos_RCP.hpp"

// Forward declarations
namespace NOX {
  class Utils;
  class GlobalData;
  class Observer;
  namespace MeritFunction {
    class Generic;
  }
  namespace Direction {
    class Generic;
  }
}

namespace NOX {
namespace Solver {

/*!
  \brief %Newton-like solver using a trust region.

Our goal is to solve: \f$ F(x) = 0, \f$ where \f$ F:\Re^n \rightarrow
\Re^n \f$.  Alternatively, we might say that we wish to solve

&nbsp;&nbsp;
\f$
\min f(x) \equiv \frac{1}{2} \|F(x)\|^2_2.
\f$

The trust region subproblem (TRSP) at iteration \f$ k\f$ is given by

&nbsp;&nbsp;
\f$
\min \; m_k(s) \equiv f_k + g_k^T d + \frac{1}{2} d^T B_k d,
\mbox{ s.t. } \|d\| \leq \Delta_k
\quad \mbox{(TRSP)}
\f$

where

<ul>
<li>\f$ f_k = f(x_k) = \frac{1}{2} \|F(x_k)\|^2_2 \f$,
<li>\f$ g_k = \nabla f(x_k) = J(x_k)^T F(x_k) \f$,
<li>\f$ B_k =  J(x_k)^T J(x_k) \approx \nabla^2 f(x_k) \f$,
<li>\f$ J(x_k)\f$ is the Jacobian of \f$ F\f$ at \f$ x_k\f$, and
<li>\f$ \Delta_k \f$ is the trust region radius.
</ul>

The "improvement ratio" for a given step \f$ s \f$ is defined as

&nbsp;&nbsp;
\f$
\rho = \displaystyle\frac{ f(x_k) - f(x_k + d) } { m_k(0) - m_k(d) }
\f$

An iteration consists of the following steps.

<ul>
<li> Compute Newton-like direction: \f$ n\f$

<li> Compute Cauchy-like direction: \f$ c\f$

<li> If this is the first iteration, initialize \f$\Delta\f$ as
     follows: If \f$\|n\|_2 < \Delta_{\min}\f$, then \f$\Delta = 2
     \Delta_{\min}\f$; else, \f$\Delta = \|n\|_2\f$.

<li> Initialize \f$\rho = -1\f$

<li> While \f$\rho < \rho_{\min}\f$ and \f$\Delta > \Delta_{\min}\f$, do the following.

     <ul>
     <li> Compute the direction \f$ d\f$ as follows:

          <ul>

          <li> If \f$\|n\|_2 < \Delta\f$, then take a Newton step by
               setting \f$ d = n\f$

      <li> Otherwise if \f$\|c\|_2 > \Delta\f$, then take a Cauchy
               step by setting \f$ d =
               \displaystyle\frac{\Delta}{\|c\|_2} c\f$

      <li> Otherwise, take a Dog Leg step by setting
               \f$ d = (1-\gamma) c + \gamma n \f$ where
           \f$
           \gamma = \displaystyle\frac
           {-c^T a + \sqrt{ (c^Ta)^2 - (c^Tc - \Delta^2) a^Ta}}{a^Ta}
           \f$
           with \f$ a = n-c\f$.

          </ul>

     <li> Set \f$ x_{\rm new} = x + d\f$ and calculate \f$ f_{\rm new}\f$

     <li> If \f$ f_{\rm new} \geq f\f$, then \f$\rho = -1\f$ Otherwise
      \f$ \rho = \displaystyle \frac {f - f_{\rm new}} {| d^T J F
      + \frac{1}{2} (J d)^T (J d)|} \f$

     </ul>

<li> Update the solution: \f$ x = x_{\rm new}\f$

<li> Update trust region:

     <ul>

     <li> If \f$\rho < \rho_{\rm s}\f$ and \f$\|n\|_2 < \Delta\f$,
          then shrink the trust region to the size of the Newton step:
          \f$\Delta = \|n\|_2\f$.

     <li> Otherwise if \f$\rho < \rho_{\rm s}\f$, then shrink the
          trust region: \f$\Delta = \max \{ \beta_{\rm s} \Delta,
          \Delta_{\min} \} \f$.

     <li> Otherwise if \f$\rho > \rho_{\rm e}\f$ and \f$\|d\|_2 =
          \Delta\f$, then expand the trust region: \f$\Delta = \min \{
          \beta_{\rm e} \Delta, \Delta_{\rm max} \} \f$.

     </ul>

</ul>

<B>Input Paramters</B>

The following parameters should be specified in the "Trust Region"
sublist based to the solver.

- "Inner Iteration Method" - Choice of trust region algorithm to use.
  Choices are:
  <ul>
  <li> "Standard Trust Region"
  <li> "Inexact Trust Region"
  </ul>

- "Direction" - Sublist of the direction parameters for the %Newton
  point, passed to the NOX::Direction::Manager constructor. If this
  sublist does not exist, it is created by default. Furthermore, if
  "Method" is not specified in this sublist, it is added with a value of
  "Newton".

- "Cauchy %Direction" - Sublist of the direction parameters for the
  Cauchy point, passed to the NOX::Direction::Manager constructor. If
  this sublist does not exist, it is created by default. Furthermore, if
  "Method" is not specified in this sublist, it is added with a value of
  "Steepest Descent" Finally, if the sub-sublist "Steepest Descent" does
  not exist, it is created and the parameter "Scaling Type" is added and
  set to "Quadratic Min Model".

- "Minimum Trust Region Radius" (\f$\Delta_{\min}\f$) - Minimum allowable trust region
  radius. Defaults to 1.0e-6.

- "Maximum Trust Region Radius" (\f$\Delta_{\max}\f$) - Minimum allowable trust region
  radius. Defaults to 1.0e+10.

- "Minimum Improvement Ratio" (\f$\rho_{\min}\f$) - Minimum improvement ratio to accept
  the step. Defaults to 1.0e-4.

- "Contraction Trigger Ratio" (\f$\rho_{\rm s}\f$) - If the improvement ratio is less than
  this value, then the trust region is contracted by the amount
  specified by the "Contraction Factor". Must be larger than "Minimum
  Improvement Ratio". Defaults to 0.1.

- "Contraction Factor" (\f$\beta_{\rm s}\f$) - See above. Defaults to 0.25.

- "Expansion Trigger Ratio" (\f$\rho_{\rm e}\f$) - If the
  improvement ratio is greater than this value, then the trust region
  is contracted by the amount specified by the "Expansion
  Factor". Defaults to 0.75.

- "Expansion Factor"  (\f$\beta_{\rm e}\f$) - See above. Defaults to 4.0.

- "Recovery Step" - Defaults to 1.0.

- "Use Ared/Pred Ratio Calculation" (boolean) - Defaults to false. If
  set to true, this option replaces the algorithm used to compute the
  improvement ratio, \f$ \rho \f$, as described above.  The improvement
  ratio is replaced by an "Ared/Pred" sufficient decrease criteria
  similar to that used in line search algorithms (see Eisenstat and
  Walker, SIAM Journal on Optimization V4 no. 2 (1994) pp 393-422):
  - \f$\rho = \frac{\|F(x) \| - \| F(x + d) \| }
                   {\| F(x) \| - \| F(x) + Jd \| } \f$

- "Use Cauchy in Newton Direction" - Boolean.  Used only by the
  "Inexact Trust Region" algorithm.  If set to true, the initial guess
  for the Newton direction computation will use the Cauchy direction as
  the initial guess.  Defaults to false.

- "Use Dogleg Segment Minimization" - Boolean.  Used only by the
  "Inexact Trust Region" algorithm.  If set to true, the \f$ \tau \f$
  parameter is minimized over the dogleg line segments instead of
  being computed at the trust regioin radius.  Used only by the
  "Inexact Trust Region" algorithm.  Defaults to false.

- "Use Counters" - Boolean.  If set to true, solver statistics will be stored.  Defaults to true.

- "Write Output Parameters" - Boolean.  If set to true, the solver statistics will be written to the relevant "Output" sublists (see Output Parameters).  Defaults to true.

- "Solver Options" - Sublist of general solver options.
   <ul>
   <li> "User Defined Pre/Post Operator" is supported.  See NOX::Parameter::PrePostOperator for more details.
   </ul>

<B>Output Paramters</B>

A sublist called "Output" will be created at the top level of the parameter list and contain the following general solver parameters:

- "Nonlinear Iterations" - Number of nonlinear iterations

- "2-Norm or Residual" - Two-norm of final residual

A sublist called "Output" will be created in the "Trust Region" sublist and contain the following trust region specific output parameters:

- "Number of Cauchy Steps" - Number of cauchy steps taken during the solve.

- "Number of Newton Steps" - Number of Newton steps taken during the solve.

- "Number of Dogleg Steps" - Number of Dogleg steps taken during the solve.

- "Number of Trust Region Inner Iterations" - Number of inner iterations
  required to adjust the trust region radius.

- "Dogleg Steps: Average Fraction of Newton Step Length" - Average value of the fraction a dogleg step took compared to the full Newton step.  The fractional value is computed as \f$ \mbox{frac} = \frac{\| d \|}{\| n\|} \f$.

- "Dogleg Steps: Average Fraction Between Cauchy and Newton Direction" -  Average value of the fraction a dogleg step took between the Cauchy and Newton directions.  This is the \f$ \gamma \f$ variable in the standard dogleg algorithm and the \f$ \tau \f$ parameter in the inexact dogleg algorithm.  A value of 0.0 is a full step in the Cauchy direction and a value of 1.0 is a full step in the Newton direction.

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class InexactTrustRegionBased : public Generic {

public:

  /*!
    \brief Constructor

    See reset() for description.
  */
  InexactTrustRegionBased(
     const Teuchos::RCP<NOX::Abstract::Group>& grp,
     const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
     const Teuchos::RCP<Teuchos::ParameterList>& params);

  //! Destructor
  virtual ~InexactTrustRegionBased();

  virtual void reset(const NOX::Abstract::Vector& initialGuess,
             const Teuchos::RCP<NOX::StatusTest::Generic>& tests);
  virtual void reset(const NOX::Abstract::Vector& initialGuess);
  virtual void reset();
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual NOX::StatusTest::StatusType getStatus() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;

  inline virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const {return solnPtr;};
  inline virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const {return oldSolnPtr;};
  inline virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const {return paramsPtr;};
  virtual Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const;

protected:

  //! "Standard" trust region implementation
  virtual NOX::StatusTest::StatusType iterateStandard();
  //! "Inexact Trust Region"
  virtual NOX::StatusTest::StatusType iterateInexact();

  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

  //! Print an error message and throw an error during parameter reads
  virtual void invalid(const std::string& param, double value) const;

  //! Print an error message and throw an error
  virtual void throwError(const std::string& method, const std::string& mesage) const;

  //! Check to see if the current step is acceptable.  If not, it reduces the trust region radius accordingly.
  NOX::StatusTest::StatusType checkStep(const NOX::Abstract::Vector& step,
                    double& radius);

  //! Computes the norm of a given vector.
  /*! Defaults to the L-2 norm but could use a user defined norm also. */
  virtual double computeNorm(const NOX::Abstract::Vector& v);

protected:

  //! Type of Trust Region algorithm to use.
  enum TrustRegionType {
    //! Basic trust region method for nonlinear systems (Nocedal and Wright?).
    Standard,
    //! Inexact Trust region WITHOUT minimization of the local linear model over the line segments.
    Inexact
  };

  //! Type of trust region algorithm to use.
  TrustRegionType method;

  //! Return types for inner iteration status test
  enum InnerIterationReturnType {
    //! Converged.
    Converged,
    //! Unconverged.
    Unconverged,
    //! Failed by hitting minimum radius bound.
    Failed
  };

  //! Pointer to the global data object.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Utils
  Teuchos::RCP<NOX::Utils> utils;

  //! Current status of the trust region inner iteration.
  InnerIterationReturnType innerIterationStatus;

  //! Current solution.
  Teuchos::RCP<NOX::Abstract::Group> solnPtr;

  //! Previous solution pointer.
  Teuchos::RCP<NOX::Abstract::Group> oldSolnPtr;

  //! Current newton direction pointer.
  Teuchos::RCP<NOX::Abstract::Vector> newtonVecPtr;

  //! Current cauchy direction pointer.
  Teuchos::RCP<NOX::Abstract::Vector> cauchyVecPtr;

  //! Extra vector used in computations
  Teuchos::RCP<NOX::Abstract::Vector> rCauchyVecPtr;

  //! Extra vector used in computations
  Teuchos::RCP<NOX::Abstract::Vector> residualVecPtr;

  //! Extra vector used in computations
  Teuchos::RCP<NOX::Abstract::Vector> aVecPtr;

  //! Extra vector used in computations
  Teuchos::RCP<NOX::Abstract::Vector> bVecPtr;

  //! Stopping test.
  Teuchos::RCP<NOX::StatusTest::Generic> testPtr;

  //! Input parameters.
  Teuchos::RCP<Teuchos::ParameterList> paramsPtr;

  //! Inexact Newton utitilities.
  NOX::Direction::Utils::InexactNewton inNewtonUtils;

  //! %Newton %Search %Direction.
  Teuchos::RCP<NOX::Direction::Generic> newtonPtr;

  //! Cauchy %Search %Direction.
  Teuchos::RCP<NOX::Direction::Generic> cauchyPtr;

  //! Radius of the trust region
  double radius;

  //! Minimum improvement ratio to accept step
  double minRatio;

  //! Minimum trust region radius
  double minRadius;

  //! Maximum trust region radius
  double maxRadius;

  //! ratio < alpha triggers contraction
  double contractTriggerRatio;

  //! ratio > beta triggers expansion
  double expandTriggerRatio;

  //! Expansion factor
  double expandFactor;

  //! Constraction factor
  double contractFactor;

  //! Take a step of this length in the Newton direction if the
  //! trust-region search fails
  double recoveryStep;

  //! Value of \f$ f \f$ at current solution
  double newF;
  //! Value of \f$ f \f$ at previous solution
  double oldF;

  //! norm(xnew - xold)
  double dx;

  //! Number of nonlinear iterations.
  int nIter;

  //! Current linear solve tolerance (inexact only).
  double eta;

  //! Linear solve tolerance used in last iteration (inexact only).
  double eta_last;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Type of check to use for status tests.  See NOX::StatusTest for more details.
  NOX::StatusTest::CheckType checkType;

  //! Enumerated list for each direction that may be required in the Trust region computation.
  enum StepType
  {
    //! Use the Newton direction
    Newton,
    //! Use the Cauchy direction
    Cauchy,
    //! Use the doglog direction
    Dogleg
  };

  //! Type of step to be taken.
  StepType stepType;

  //! Stores merit function supplied by global data.
  Teuchos::RCP<NOX::MeritFunction::Generic> meritFuncPtr;

  //! If set to true, the initial guess for the Newton direction computation will use the Cauchy direction as the initial guess.
  bool useCauchyInNewtonDirection;

  //! If set to true, statistics/counters will be output to the output list
  bool writeOutputParamsToList;

  //! If set to true, counters will be stored by the solver.
  bool useCounters;

  //! Counters for the algorithm
  NOX::SolverStats::TrustRegionStats* counters;

  //! If set to true, the minimum improvement ratio condition uses an Ared/Pred approach.
  bool useAredPredRatio;

  //! If set to true, the \f$ \tau \f$ parameter is minimized over the dogleg line segments instead of being computed at the trust regioin radius.
  bool useDoglegMinimization;

  //! Pointer to a user defined NOX::Observer object.
  Teuchos::RCP<NOX::Observer> observer;

};
} // namespace Solver
} // namespace NOX

#endif
