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

#ifndef NOX_DIRECTION_NEWTON_H
#define NOX_DIRECTION_NEWTON_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Parameter_List.H"	   // class data element 

// Forward declaration
namespace NOX {
class Utils;
}

namespace NOX {
namespace Direction {

//! %Newton direction computation
/*! 
  Computes the Newton direction by solving
  \f[ Js = -F. \f]
  
  To use this direction, specify that the "Method" is "Newton" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  If the underlying linear solver for this comupation is an iterative method, 
  this direction computation becomes an inexact %Newton method.  The 
  tolerance used in terminating the linear solve (called a \b Forcing \b Term) 
  may be constant or it can be adjusted for each successive linear solve.
  Adjustable Forcing Terms are implemented as described in [1] and
  originally introduced in [2]. In this case, we choose a convergence
  tolerance for the linear solve. In other words, at nonlinear iteration 
  \f$k\f$ we choose \f$\eta_k \in \left[ \eta_{\min}, \eta_{\max} \right] \f$
  and solve for \f$s_k\f$ such that

  \f[
  \frac{\|F_k + J_k s_k \|}{\|F_k\|} \leq \eta_k
  \f]

  The \b Type \b 1 method sets the forcing term as follows:

  \f[
  \eta_k = \left\vert \frac{\| F_k \| - \|F_{k-1} + J_{k-1} s_{k-1} \| }
  {\|F_{k-1}\|} \right\vert
  \f]

  With the following safeguards imposed:
  \f[
  \max\{\eta_{k-1}^{\frac{1 + \sqrt{5}}{2}}, \eta_{\min} \} \leq \eta_k \leq \eta_{\max}
  \f]

  The \b Type \b 2 method sets the forcing term as follows:

  \f[
  \eta_k = \gamma \left( \frac{\|F_k\|}{\|F_{k-1}\|} \right)^\alpha
  \f]

  With the following safeguards imposed:
  \f[
  \max\{\gamma \eta_{k-1}^{\alpha}, \eta_{\min} \} \leq \eta_k \leq \eta_{\max}
  \f]

  The following options may be specified in the "Newton" sublist of
  the "Direction" sublist of the solver parameters.

  - "Rescue Bad Newton Solve" (Boolean) - Says whether or not we should
    try to use the solution from a bad linear solve. If set to true, 
    we will accept the solution if it produces any reduction whatsoever
    in the linear residual. Defaults to true.

  - "Forcing Term Method" - Method to compute the forcing term, i.e.,
    the tolerance for the linear solver. Choices are "Constant", 
    "Type 1", and "Type 2".  Defaults to "Constant". 

  - "Forcing Term Initial Tolerance" - Initial linear
    solver tolerance (\f$\eta_0\f$). Defaults to 1.0e-1.

  - "Forcing Term Minimum Tolerance" - Minimum acceptable linear
    solver tolerance (\f$\eta_{\min}\f$). Defaults to 1.0e-6.

  - "Forcing Term Maximum Tolerance" = Maximum acceptable linear
    solver tolerance (\f$\eta_{\max}\f$). Default to 0.01.

  - "Forcing Term Alpha" - Alpha (\f$\alpha\f$) used for the "Type 2"
    forcing term calcuation. Defaults to 1.5.

  - "Forcing Term Gamma" - Gamma (\f$\gamma\f$) used for the "Type 2"
    forcing term calcuation. Defaults to 0.9.

  <B> Important Note: Using inexact %Newton methods with adjustable forcing terms </B> <br>
The performance of adjustable forcing terms can be affected if the linear solver solves a linear system that is scaled.  The forcing term algorithm assumes that the linear solver achieves the following tolerance:

  \f[
  \frac{\|F_k + J_k s_k \|}{\|F_k\|} \leq \eta_k
  \f]

However, if we are solving a scaled linear system we are actually solving to the following criteria:

  \f[
  \frac{\|D_k^{-1}F_k + D_k^{-1}J_k s_k \|}{\|D_k^{-1}F_k\|} \leq \eta_k
  \f]

Where \f$ D_k \f$ is the scaling matrix (typically a diagonal row sum matrix is used).  This is a different tolerance than the one expected in the adjustable forcing term algorithm.  To account for this, one would have to have the NOX::Direction::Newton compute \f$ D \f$-norms in the direction computations.  This breaks the abstract design (and bloats the NOX::Abstract::Group) since now NOX::Direction::Newton would have to know details about the linear solver and then access scaling objects to compute these special norms.  Alternatively, we have chosen to change the convergence test in the linear solver to achieve the same goal.  The linear solver must now account for scaling in the termination criteria.  In %NOX, if an adjustable forcing term is used in conjunction with scaling of the linear system by a scaling matrix \f$ D \f$, the following convergence criteria should instead be used to terminate the linear solve:

  \f[
  \frac{\| D\left[D^{-1}F_k + D^{-1}J_k s_k \right] \|}{\|F_k\|} \leq \eta_k
  \f]

NOX::Direction::Newton will set a boolean parameter called "Using Adjustable Forcing Term" in the "%Linear %Solver" sublist of the "%Newton" direction sublist.  This parameter will be set to \b true if the "Forcing Term Method" is set to "Type 1" or "Type 2".  It is set to \b false if the "Forcing Term Method" is set to "Constant".  A linear solver can check this parameter and thus adjust it's convergence criteria accordingly.

  \b References

  <ol>

  <li>Michael Pernice and Homer F. Walker, NITSOL: A %Newton Iterative
  %Solver for Nonlinear Systems, SISC 19(Jan 1998):302-318.

  <li> S. C. Eisenstat and H. F. Walker, Globally convergent inexact
  %Newton methods, SINUM 19(1982):400-408

  </ol>
 */  
class Newton : public Generic {

public:

  //! Constructor
  Newton(const NOX::Utils& u, NOX::Parameter::List& params);

  //! Destructor
  virtual ~Newton();

  // derived
  virtual bool reset(NOX::Parameter::List& params);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::Generic& solver);

protected:

  //! Called each iteration to reset the forcing term (ie, the
  //! convergence tolerance for the linear solver).
  virtual bool resetForcingTerm(const NOX::Abstract::Group& soln, 
				const NOX::Abstract::Group& oldSoln, 
				int niter, 
				const NOX::Parameter::List& solverParams);


private:

  //! Printing Utils
  const NOX::Utils& utils;

  //! Return true if we can recover from a bad Newton solve
  /*!
    In other words, sometimes we'll use the calculated Newton
    direction even if it doesn't satisfy the convergence tolerance.

    If doRescue is true,
    \f[
    \|J s - (-F) \| / \|F\| < 1,
    \f]
    and we don't encounter any other errors, return true.
  */
  bool rescueBadNewtonSolve(const NOX::Abstract::Group& soln) const;

  //! Print an error message
  void throwError(const string& functionName, const string& errorMsg);


private:

  //! "Direction" sublist with parameters for the direction vector
  /*! 
    \note This is pointer rather than a reference to allow for the
    reset function. 
   */
  NOX::Parameter::List* paramsPtr;

  //! Determined based on "Rescue Bad Newton Solve". Used in rescueBadNewtonSolve
  bool doRescue;

  //! Vector containing the predicted RHS
  NOX::Abstract::Vector* predRhs;

  //! Vector containing the computed direction
  NOX::Abstract::Vector* stepDir;

};
} // namespace Direction
} // namespace NOX
#endif
