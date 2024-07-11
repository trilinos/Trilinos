// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_DIRECTION_NEWTON_H
#define NOX_DIRECTION_NEWTON_H

#include "NOX_Direction_Generic.H" // base class
#include "Teuchos_ParameterList.hpp"       // class data element
#include "Teuchos_RCP.hpp"

// Forward declaration
namespace NOX {
  using std::string;
  class Utils;
  namespace Parameter {
    class UserNorm;
  }
}

namespace NOX {
namespace Direction {

//! %Newton direction computation
/*!
  Computes the %Newton direction by solving the %Newton system.
  \f[
  Jd = -F
  \f]

  Here \f$J\f$ is the n x n Jacobian matrix at the current
  iterate, \f$F\f$ is the n-vector representing the nonlinear
  function at the current iterate, and \f$d\f$ is the n-vector
  that we are solving for.

  If we use an iterative linear solver for the %Newton system, then
  this is called an <b>inexact %Newton method</b>.  The tolerance used
  to terminate the linear solve is called the <b>forcing term</b>.
  The forcing term may be constant, or it may be <b>adjustable</b>.
  In either case, at iteration \f$k\f$ we require,
  \f[
  \frac{\|J_k d_k - (-F_k)\|}{\|F_k\|} \leq \eta_k.
  \f]
  Here \f$\eta_k\f$ is the forcing term for iteration \f$k\f$.

  \note This solution tolerance is to be enforced by the user's
  implementation of NOX::Abstract::Group::computeNewton; it is passed
  in as the "Tolerance" in the parameter list for that function.

  Adjustable forcing terms were introduced by Eisenstat and Walker
  (1982); here they are implemented as described in Pernice and Walker
  (1998). We have two choices for adjustable forcing terms:

  <ul>
  <li><b>Type 1</b>

  \f[
  \eta_k = \left\vert \frac{\| F_k \| - \|J_{k-1} d_{k-1} - (-F_{k-1}) \| }
  {\|F_{k-1}\|} \right\vert
  \f]

  With the following safeguards imposed:
  \f[
  \max\{\eta_{k-1}^{\frac{1 + \sqrt{5}}{2}}, \eta_{\min} \} \leq \eta_k \leq \eta_{\max}
  \f]

  <li><b>Type 2</b>

  \f[
  \eta_k = \gamma \left( \frac{\|F_k\|}{\|F_{k-1}\|} \right)^\alpha
  \f]

  With the following safeguards imposed:
  \f[
  \max\{\gamma \eta_{k-1}^{\alpha}, \eta_{\min} \} \leq \eta_k \leq \eta_{\max}
  \f]
  </ul>


<B>Parameters</B>

"Direction":
<ul>
<li> "Method" = "Newton" [required]
</ul>

"Direction"/"Newton":

<ul>
<li>"Forcing Term Method" - Method to compute the forcing term, i.e.,
    the tolerance for the linear solver. Choices are:
    <ul>
    <li>"Constant" [default]
    <li>"Type 1"
    <li>"Type 2"
    </ul>

<li>"Forcing Term Initial Tolerance" - \f$\eta_0\f$ (initial linear
    solver tolerance). Defaults to 0.1.

<li>"Forcing Term Minimum Tolerance" - \f$\eta_{\min}\f$. Defaults to 1.0e-6.

<li>"Forcing Term Maximum Tolerance" - \f$\eta_{\max}\f$. Defaults to 0.01.

<li>"Forcing Term Alpha" - \f$\alpha\f$ (used only by "Type 2"). Defaults to 1.5.

<li>"Forcing Term Gamma" - \f$\gamma\f$ (used only by "Type 2"). Defaults to 0.9.

<li>"Rescue Bad %Newton Solve" (Boolean) - If set to true, we will use
the computed direction even if the linear solve does not achieve the
tolerance specified by the forcing term. Defaults to true.

</ul>

"Direction"/"Newton"/"Linear Solver":

<ul>

<li> "Tolerance" - Tolerance for the linear solve. This may be
adjusted automatically by the forcing calculation. Defaults to
1.0e-10. <em>Will be adjusted automatically by NOX if the "Forcing
Term Method" is "Type 1" or "Type 2".</em>

</ul>

\note
When using a forcing term, it's critically important the the residual
of the original system is used in the comparison. This can be an issue
if scaling or left preconditioning is applied to the linear system.

<b>References</b>

  <ul>

  <li>Michael Pernice and Homer F. Walker, <em>NITSOL: A %Newton Iterative
  %Solver for Nonlinear Systems</em>, SISC 19(Jan 1998):302-318.

  <li> S. C. Eisenstat and H. F. Walker, <em>Globally convergent inexact
  %Newton methods</em>, SINUM 19(1982):400-408

  </ul>
*/


class Newton : public Generic {

public:

  //! Constructor
  Newton(const Teuchos::RCP<NOX::GlobalData>& gd,
     Teuchos::ParameterList& params);

  //! Destructor
  virtual ~Newton();

  // derived
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& params);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::LineSearchBased& solver);

protected:

  //! Called each iteration to reset the forcing term (ie, the
  //! convergence tolerance for the linear solver).
  virtual bool resetForcingTerm(const NOX::Abstract::Group& soln,
                const NOX::Abstract::Group& oldSoln,
                int niter,
                const NOX::Solver::Generic& solver);


private:

  //! Print an error message
  void throwError(const std::string& functionName, const std::string& errorMsg);


private:

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Printing Utilities
  Teuchos::RCP<NOX::Utils> utils;

  //! "Direction" sublist with parameters for the direction vector
  /*!
    \note This is pointer rather than a reference to allow for the
    reset function.
   */
  Teuchos::ParameterList* paramsPtr;

  //! Determined based on "Rescue Bad %Newton Solve".
  bool doRescue;

  //! Vector containing the predicted RHS
  Teuchos::RCP<NOX::Abstract::Vector> predRhs;

  //! Vector containing the computed direction
  Teuchos::RCP<NOX::Abstract::Vector> stepDir;

  //! Flag used to determine if adjustable forcing term is being used.
  bool useAdjustableForcingTerm;

  //! Current linear solve tolerance.
  double eta_k;

  //! Minimum linear solve tolerance.
  double eta_min;

  //! Maximum linear solve tolerance.
  double eta_max;

  //! Initial linear solve tolerance.
  double eta_initial;

  //! Name of the method used for calculating the inexact forcing term.
  std::string method;

  //! Parameter used for Type 2 forcing term calculation.
  double alpha;

  //! Parameter used for Type 2 forcing term calculation.
  double gamma;
};
} // namespace Direction
} // namespace NOX
#endif
