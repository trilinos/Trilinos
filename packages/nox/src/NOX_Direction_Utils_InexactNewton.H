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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_DIRECTION_UTILS_INEXACTNEWTON_H
#define NOX_DIRECTION_UTILS_INEXACTNEWTON_H

#include "NOX_Common.H"  // for std::string
#include "Teuchos_RCP.hpp"

// Forward declaration
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class Utils;
  class GlobalData;
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Solver {
    class Generic;
  }
}

namespace NOX {
namespace Direction {
namespace Utils {

//! %Inexact Newton Utilities
/*! 

  If we use an iterative linear solver for a %Newton-based  solve, then
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

<li>"Forcing Term User Defined Norm" (NOX::Parameter::UserNorm derived
    object) - If using a Type 1 or Type 2 adjustable forcing term, the
    norms used to calculate \f$ \eta_k \f$ should be based on the same
    norm used in the convergence test of the linear solver.
    Essentially this means that the norm must account for LEFT scaling
    of any kind that is applied to the linear system during a solve.
    If set, the computation of \f$ \eta_k \f$ will be done using a
    used defined function that is passed in through a
    NOX::Parameter::Arbitrary derived object. It will take the
    arbitrary object and cast it to a NOX::Parameter::UserNorm object
    and use it for norm computations.  If this parameter is not set,
    this method uses the L-2 Norm for any norm computations of \f$
    \eta_k \f$.

<li>"Set Tolerance in Parameter List" When calling
computeForcingTerm() the value of the forcing term will be set in the
parmaeter list pointed to by the InexactNewton object.  It will be set
under paramsPtr->sublist(<directionMethod>).sublist("Linear
Solver").set("Tolerance", eta_k). Defaults to true.

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


class InexactNewton {

public:

  //! Constructor
  InexactNewton(const Teuchos::RCP<NOX::GlobalData>& gd,
		Teuchos::ParameterList& directionSublist);

  //! Destructor
  virtual ~InexactNewton();

  //! Reset the utilities
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd, 
		     Teuchos::ParameterList& directionSublist);

  //! Called each iteration to reset the forcing term (ie, the
  //! convergence tolerance for the linear solver).
  /*!  if the user supplied eta_last then it will use this value for
    eta_km1 instead of looking for it in the "Linear Solver" sublist.
  */
  virtual double computeForcingTerm(const NOX::Abstract::Group& soln, 
				    const NOX::Abstract::Group& oldSoln, 
				    int niter, 
				    const NOX::Solver::Generic& solver,
				    double eta_last = -1.0);

private:

  //! Print an error message
  void throwError(const std::string& functionName, const std::string& errorMsg);

private:

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Printing Utilities
  Teuchos::RCP<NOX::Utils> printing;

  //! Type of forcing term calculations available
  enum ForcingTermType {
    //! Constant
    Constant,
    //! Type 1
    Type1,
    //! type 2
    Type2
  };

  //! Type of forcing term calculation to use
  ForcingTermType forcingTermMethod;

  //! "Inexact Newton" sublist with parameters for the direction vector
  /*! 
    \note This is a pointer rather than a reference to allow for the
    reset function. 
   */
  Teuchos::ParameterList* paramsPtr;

  //! Vector containing the predicted RHS
  Teuchos::RCP<NOX::Abstract::Vector> predRhs;

  //! Vector containing the step direction
  Teuchos::RCP<NOX::Abstract::Vector> stepDir;

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
  
  //! Name of the direction method that forcing terms will be used by.
  std::string directionMethod;
  
  //! Parameter used for Type 2 forcing term calculation.
  double alpha;

  //! Parameter used for Type 2 forcing term calculation.
  double gamma;

  //! Determines whether to set the "Tolerance" in the parameter list.
  bool setTolerance;
};
} // namespace Utils
} // namespace Direction
} // namespace NOX
#endif
