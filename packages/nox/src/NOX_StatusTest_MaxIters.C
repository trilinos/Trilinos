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

NOX::StatusTest::MaxIters::MaxIters(int maxIterations) :
  maxiters(maxIterations),
  status(Unevaluated)
{
  if (maxiters < 1) 
  {
    cout << "NOX::StatusTest::MaxIters - must choose a number greater than zero" << endl;
    throw "NOX Error";
  }
}

NOX::StatusTest::MaxIters::~MaxIters()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::MaxIters::checkStatus(const Solver::Generic& problem)
{
  return checkStatus(problem, NOX::StatusTest::Minimal);
}

NOX::StatusTest::StatusType NOX::StatusTest::MaxIters::checkStatus(const Solver::Generic& problem, 
								   NOX::StatusTest::CheckType checkType)
{
  switch (checkType)
  {
  case NOX::StatusTest::Complete:
  case NOX::StatusTest::Minimal:
    niters = problem.getNumIterations();
    status = (niters >= maxiters) ? Failed : Unconverged;
    break;

  case NOX::StatusTest::None:
  default:
    niters = -1;
    status = Unevaluated;
    break;
  }

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::MaxIters::getStatus() const
{
  return status;
}

ostream& NOX::StatusTest::MaxIters::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Number of Iterations = " << niters << " < " << maxiters;
  stream << endl;
 return stream;
}

int NOX::StatusTest::MaxIters::getMaxIters() const
{
  return maxiters;
}

int NOX::StatusTest::MaxIters::getNumIters() const
{
  return niters;
}
