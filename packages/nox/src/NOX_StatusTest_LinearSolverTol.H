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

#ifndef NOX_STATUSTEST_LINEARSOLVERTOLERANCE_H
#define NOX_STATUSTEST_LINEARSOLVERTOLERANCE_H

#include "NOX_StatusTest_Generic.H"	// base class

namespace NOX {

namespace StatusTest {

//! Convergence test that forces the linear solver tolerance on a converged iteration to be below a value specified in the constructor (defaults to 0.5).
/*! 
  There are instances where the linear solver fails to make any progress.  This occurs when (1) the adjustable convergence criteria is used and the linear solver tolerance computed is too large or (2) the linear solver fails to achieve the requested convergence yet returns its best approximation to the solution.  In these cases, the change in the solution vector between iterates becomes very small and if we are using, NormUpdate and/or NormWRMS as the convergence test, the solution can be prematurely flagged as converged.  To prevent this from occurring, we use this StatusTest to make sure the linear solver is making sufficient progress.

NOTE: This status test requires the linear solver to set the final scaled tolerance in the "Output" sublist of the "Linear Solver" sublist of the "Direction" sublist: 

parameterList.sublist("Direction").sublist("Linear Solver").sublist("Output").setParameter("Scaled Residual", value);

*/
class LinearSolverTol : public Generic {

public:

  //! Constructor. 
  LinearSolverTol(double tolerance = 0.5);

  //! Destructor.
  virtual ~LinearSolverTol();

  virtual StatusType checkStatus(const Solver::Generic& problem);

  virtual StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  //! %Status
  StatusType status;

  //! Requested linear solver tolerance.
  double specifiedTolerance;

  //! Actual linear solver tolerance.
  double tolerance;

};

} // namespace Status
} // namespace NOX

#endif
