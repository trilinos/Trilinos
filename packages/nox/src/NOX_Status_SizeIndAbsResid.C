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

#include "NOX_Status_SizeIndAbsResid.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX::Status;

SizeIndAbsResid::SizeIndAbsResid(double tolerance)
{
  tol = tolerance;
  status = Unconverged;
}

SizeIndAbsResid::~SizeIndAbsResid()
{
}

StatusType SizeIndAbsResid::operator()(const Solver::Generic& problem)
{
  status = Unconverged;
  const Abstract::Group& tmp = problem.getSolutionGroup();
  double normrhs = tmp.getNormRHS();
  normrhs /= sqrt( (double) tmp.getRHS().length() );

  //cout << "Size Ind Residual 2-Norm = " << normrhs << endl; 

  if (normrhs < tol)
    status = Converged;
  return status;
}

ostream& SizeIndAbsResid::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Size Independent Residual 2-Norm < " << Utils::sci(tol);
  stream << endl;
  return stream;
}
