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

#include "NOX_Common.H"
#include "NOX_Direction_SteepestDescent.H" // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"

using namespace NOX;
using namespace NOX::Direction;

SteepestDescent::SteepestDescent(const Parameter::List& params) 
{
  reset(params);
}

SteepestDescent::~SteepestDescent()
{
  
}

bool SteepestDescent::reset(const Parameter::List& params)
{
  return true;
}

bool SteepestDescent::operator()(Abstract::Vector& dir, 
				 Abstract::Group& soln, 
				 const Solver::Generic& solver) 
{
  // Compute RHS at current solution
  bool ok = soln.computeRHS();

  if (!ok) {
    cerr << "NOX::Direction::SteepestDescent::operator() - Unable to compute RHS." << endl;
    throw "NOX Error";
  }
  
  // Compute Jacobian at current solution.
  ok = soln.computeJacobian();

  if (!ok) {
    cerr << "NOX::Direction::SteepestDescent::operator() - Unable to compute Jacobian." << endl;
    throw "NOX Error";
  }
  
  // Compute the gradient
  ok = soln.computeGrad();

  if (!ok) {
    cerr << "NOX::Direction::SteepestDescent::operator() - Unable to compute gradient." << endl;
    throw "NOX Error";
  }
  
  // Get the gradient direction.
  dir = soln.getGrad();

  // Compute the steepest descent direction
  double norm = dir.norm();
  dir.scale(-1.0/norm);

  return true;
}


