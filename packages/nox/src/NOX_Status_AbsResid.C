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

#include "NOX_Status_AbsResid.H"

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"

using namespace NOX::Status;

AbsResid::AbsResid(double tolerance)
{
  tol = tolerance;
  status = Unconverged;
}

AbsResid::~AbsResid()
{
}

StatusType AbsResid::operator()(const Solver::Generic& problem)
{
  const Abstract::Group& tmp = problem.getSolutionGroup();
  double normrhs = tmp.getNormRHS();
  if (normrhs < tol)
    status = Converged;
  return status;
}

ostream& AbsResid::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "Absolute Residual Norm with Tolerance < " << tol;
  stream << " : " << ((status == Unconverged) ? "Unconverged" : "CONVERGED!");
  stream << endl;
  return stream;
}
