// $Id$ 
// $Source$ 

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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_DIRECTION_MODIFIEDNEWTON_H
#define NOX_DIRECTION_MODIFIEDNEWTON_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Parameter_List.H"	   // class data element 


namespace NOX {
  class Utils;                     // Forward declaration

namespace Direction {

//! %Modified Newton direction computation
/*! 
  Computes the Newton direction by solving
  \f[ Js = -f. \f]
  where \f[ J \f] may be an old Jacobian. 
  
  To use this direction, specify that the "Method" is "Modified-Newton" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).  The default maximum age of the Jacobian is 10, and
  can be set using the "Max Age of Jacobian" parameter in the "Modified-Newton"
  sublist.  
  Max Age of Jacobian  = 0 => Secant method (never update Jacobian)
                       = 1 => Newton's method (always update Jacobian)
                      >= 2 => Modified Newton's method

  \b References

  <ol>

  <li>C. T. Kelley, Iterative methods for linear and nonlinear equations, volume 16 of Frontiers in Applied Mathematics, SIAM, Philadelphia PA, 1995.

  </ol>
 */  
class ModifiedNewton : public Generic {

public:

  //! Constructor
  ModifiedNewton(const NOX::Utils& u, NOX::Parameter::List& params);

  //! Destructor
  virtual ~ModifiedNewton();

  // derived
  virtual bool reset(NOX::Parameter::List& params);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::Generic& solver);

protected:

private:
 
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

  //! Printing Utilities
  const NOX::Utils& utils;

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

  //! Group containing the old Jacobian
  Abstract::Group* oldJacobianGrpPtr;

  //! Integer counting how old the Jacobian is
  int ageOfJacobian;

  //! Integer containing the maximum age of the Jacobian
  int maxAgeOfJacobian;

};
} // namespace Direction
} // namespace NOX
#endif // NOX_DIRECTION_MODIFIEDNEWTON_H
#endif // WITH_PRERELEASE
