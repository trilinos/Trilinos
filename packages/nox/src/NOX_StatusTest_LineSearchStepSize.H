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

#ifndef NOX_STATUSTEST_LINESEARCHSTEPSIZE_H
#define NOX_STATUSTEST_LINESEARCHSTEPSIZE_H

#include "NOX_StatusTest_Generic.H"	// base class

namespace NOX {

namespace StatusTest {

//! Convergence test that forces the step size computed in a line search to be greater than a specified value to be flagged as converged (defaults to 1.0).
/*! 
  Line search algorithms can compute very small step sizes.  This leads to very small changes in the solution vector, \f$\Delta x \f$.   If all other convergence criteria is based on this \f$\Delta x\f$ (if we are using, for example, NormUpdate and/or NormWRMS) then the code could prematurely accept an unconverged solution as converged since the norm of the change in x is small.  To prevent this from occurring, we use this StatusTest to make sure the backtracked step length is not small.  Essentially we want to be taking full steps or close to full steps at a truly converged solution.

NOTE: This test is only valid if using the NOX::Solver::LineSearchBased solver.  Otherwise it returns a converged status.

*/
class LineSearchStepSize : public Generic {

public:

  //! Constructor. 
  LineSearchStepSize(double specifiedStepSize = 1.0);

  //! Destructor.
  virtual ~LineSearchStepSize();

  virtual StatusType checkStatus(const Solver::Generic& problem);

  virtual StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  //! %Status
  StatusType status;

  //! Minimum step size that is allowed for convergence.
  double minStepSize;

  //! Actual step size computed in the line search.
  double computedStepSize;

};

} // namespace Status
} // namespace NOX

#endif
