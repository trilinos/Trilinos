#ifdef WITH_PRERELEASE

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

#ifndef NOX_DIRECTION_MODIFIEDNEWTON_H
#define NOX_DIRECTION_MODIFIEDNEWTON_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Parameter_List.H"	   // class data element 

namespace NOX {
namespace Direction {

//! %Modified Newton direction computation
/*! 
  Computes the Newton direction by solving
  \f[ Js = -f. \f]
  
  To use this direction, specify that the "Method" is "Modified Newton" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  We may optionally use a \b Forcing \b Term, as described in [1] and
  originally introduced in [2]. In this case, we choose a convergence
  tolerance for the linear solve. In other words, at iteration \f$k\f$
  we choose \f$\eta_k \in \left[ \eta_{\min}, \eta_{\max} \right] \f$
  and solve for \f$s_k\f$ such that

  \f[
  \frac{\|F_k + J_k^T s_k \|}{\|F_k\|} \leq \eta_k
  \f]

  The \b Type \b 1 methods sets the forcing term as follows:

  \f[
  \eta_k = \left\vert \frac{\| F_k \| - \|F_{k-1} + J_{k-1} s_{k-1} \| }
  {\|F_{k-1}\|} \right\vert
  \f]

  We impose the following safeguards.
  \f[
  \max\{\eta_{k-1}^{\frac{1 + \sqrt{5}}{2}}, \eta_{\min} \} \leq \eta_k \leq \eta_{\max}
  \f]

  The \b Type \b 2 methods sets the forcing term as follows:

  \f[
  \eta_k = \gamma \left( \frac{\|F_k\|}{\|F_{k-1}\|} \right)^\alpha
  \f]

  We impose the following safeguards.
  \f[
  \max\{\gamma \eta_{k-1}^{\alpha}, \eta_{\min} \} \leq \eta_k \leq \eta_{\max}
  \f]

  The following options may be specified in the "Direction" sublist of
  the solver parameters.

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


  \b References

  <ol>

  <li>Michael Pernice and Homer F. Walker, NITSOL: A %Newton Iterative
  %Solver for Nonlinear Systems, SISC 19(Jan 1998):302-318.

  <li> S. C. Eisenstat and H. F. Walker, Globally convergent inexact
  %Newton methods, SINUM 19(1982):400-408

  </ol>
 */  
class ModifiedNewton : public Generic {

public:

  //! Constructor
  ModifiedNewton(Parameter::List& params);

  //! Destructor
  virtual ~ModifiedNewton();

  // derived
  virtual bool reset(Parameter::List& params);

  // derived
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
		       const Solver::Generic& solver);

protected:

  //! Called each iteration to reset the forcing term (ie, the
  //! convergence tolerance for the linear solver).
  virtual bool resetForcingTerm(const Abstract::Group& soln, const Abstract::Group& oldsoln, int niter, const Parameter::List& solverParams);


private:

  //! "Direction" sublist with parameters for the direction vector
  /*! 
    \note This is pointer rather than a reference to allow for the
    reset function. 
   */
  Parameter::List* paramsPtr;

  //! Vector containing the predicted RHS
  Abstract::Vector* predrhs;

  //! Vector containing the computed direction
  Abstract::Vector* stepdir;

  //! Group containing the old Jacobian
  Abstract::Group* oldJacobian;

  //! Integer counting how old the Jacobian is
  int ageOfJacobian;

};
} // namespace Direction
} // namespace NOX
#endif // NOX_DIRECTION_MODIFIEDNEWTON_H
#endif // WITH_PRERELEASE
