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

#include "NOX_StatusTest_MaxIters.H" // class definition
#include "NOX_Common.H"
#include "NOX_Solver_Generic.H"

using namespace NOX::StatusTest;

MaxIters::MaxIters(int maxiterations)
{
  if (maxiterations < 1) {
    cout << "Error: Must choose maxiterations > 1 in NOX::StatusTest::MaxIters" << endl;
    throw;
  }
    
  maxiters = maxiterations;
  status = Unconverged;
}

MaxIters::~MaxIters()
{
}

StatusType MaxIters::checkStatus(const Solver::Generic& problem)
{
  status = Unconverged;
  int niters = problem.getNumIterations();
  if (niters >= maxiters)
    status = Failed;
  return status;
}

StatusType MaxIters::getStatus() const
{
  return status;
}

ostream& MaxIters::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Number of Iterations < " << maxiters;
  stream << endl;
 return stream;
}
