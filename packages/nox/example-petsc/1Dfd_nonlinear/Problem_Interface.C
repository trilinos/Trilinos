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
                                                                                
// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "FiniteDifference.H"
using namespace std;

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteDifference& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Vec& x, Vec& RHS)
{
  return problem.evaluate(RHS_ONLY, &x, &RHS, NULL);
}

bool Problem_Interface::computeJacobian(const Vec& x, Mat& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL, &Jac);
}

bool Problem_Interface::computePreconditioner(Mat& M)
{
  cout << "ERROR: Problem_Interface::computePreconditioner() - Use Explicit Jaciban only for this test problem!" << endl;
  throw;
}
bool Problem_Interface::preconditionVector(Vec& y)
{
  cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!" << endl;
  throw;
}
//-----------------------------------------------------------------------------

