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

#include "NOX_StatusTest_LinearSolverTol.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::StatusTest;

LinearSolverTol::LinearSolverTol(double t = 0.5) :
  status(Unconverged),
  specifiedTolerance(t),
  tolerance(1.0)
{
}

LinearSolverTol::~LinearSolverTol()
{
}

StatusType LinearSolverTol::checkStatus(const Solver::Generic& problem)
{
  status = Unconverged;
  bool outputListExists = false;
  const NOX::Parameter::List& p = problem.getParameterList();

  // Make sure the output parameter list exists
  // If so, get the tolerance from it
  if (p.isParameterSublist("Direction")) {
    if (p.sublist("Direction").isParameterSublist("Linear Solver")) {
      if (p.sublist("Direction").sublist("Linear Solver")
	  .isParameterSublist("Output")) {
      
	outputListExists = true;

	tolerance = problem.getParameterList().sublist("Direction").sublist("Linear Solver").sublist("Output").getParameter("Scaled Residual", -1.0);

	if (tolerance <= specifiedTolerance)
	  status = Converged;

      }
    }
  }

  if (!outputListExists)
    return Converged;

  return status;
}

StatusType LinearSolverTol::getStatus() const
{
  return status;
}

ostream& LinearSolverTol::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Linear Solver Tolerance = " << tolerance;
  stream << " =< " << specifiedTolerance;
  stream << endl;

  return stream;
}
