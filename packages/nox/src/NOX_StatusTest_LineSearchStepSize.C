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

#include "NOX_StatusTest_LineSearchStepSize.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::StatusTest;

LineSearchStepSize::LineSearchStepSize(double s) :
  status(Unconverged),
  minStepSize(s),
  computedStepSize(0.0)
{
}

LineSearchStepSize::~LineSearchStepSize()
{
}

StatusType LineSearchStepSize::checkStatus(const Solver::Generic& problem)
{
  status = Unconverged;

  // Determine if the Generic solver is a LineSearchBased solver
  // If it is not then return a "Converged" status
  const Solver::Generic* test = 0;
  test = dynamic_cast<const Solver::LineSearchBased*>(&problem);
  if (test == 0) {
    if (Utils::doPrint(Utils::Warning)) {
      cout << "Warning: NOX::StatusTest::LineSearchStepSize::checkStatus() - "
	   << "LineSearchStepSize is only valid with the "
	   << "\"Line Search Based\" solver" << endl;
    }
  
    return Converged; 
  }

  computedStepSize = (dynamic_cast<const Solver::LineSearchBased*>(&problem))->getStepSize();

  if (computedStepSize >= minStepSize)
    status = Converged;

  return status;
}

StatusType LineSearchStepSize::getStatus() const
{
  return status;
}

ostream& LineSearchStepSize::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Line Search Step Size = " << computedStepSize;
  stream << " >= " << minStepSize;
  stream << endl;

  return stream;
}
