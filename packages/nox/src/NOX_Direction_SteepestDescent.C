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
#include "NOX_Parameter_List.H"

using namespace NOX;
using namespace NOX::Direction;

SteepestDescent::SteepestDescent(Parameter::List& params) :
  tmpVecPtr(NULL)
{
  reset(params);
}

SteepestDescent::~SteepestDescent()
{
  delete tmpVecPtr;
}

bool SteepestDescent::reset(Parameter::List& params)
{
  const string tmp = params.getParameter("Scaling Type", "2-Norm");
  if (tmp == "2-Norm")
    scaleType = SteepestDescent::TwoNorm;
  else if (tmp == "F 2-Norm")
    scaleType = SteepestDescent::FunctionTwoNorm;
  else if (tmp == "Quadratic Model Min")
    scaleType = SteepestDescent::QuadMin;
  else if (tmp == "None")
    scaleType = SteepestDescent::None;
  else {
    cout << "NOX::Direction::SteepestDescent::reset - Invalid choice \""
	  << tmp << "\" for \"Scaling Type\"" << endl;
    throw "NOX Error";
  }

 return true;
}

bool SteepestDescent::compute(Abstract::Vector& dir, 
				 Abstract::Group& soln, 
				 const Solver::Generic& solver) 
{
  // Compute F at current solution
  bool ok = soln.computeF();

  if (!ok) {
    cerr << "NOX::Direction::SteepestDescent::compute - Unable to compute F." << endl;
    throw "NOX Error";
  }
  
  // Compute Jacobian at current solution.
  ok = soln.computeJacobian();

  if (!ok) {
    cerr << "NOX::Direction::SteepestDescent::compute - Unable to compute Jacobian." << endl;
    return false;
  }
  
  // Compute the gradient
  ok = soln.computeGradient();

  if (!ok) {
    cerr << "NOX::Direction::SteepestDescent::compute - Unable to compute gradient." << endl;
    return false;
  }
  
  // Get the gradient direction.
  dir = soln.getGradient();

  // Scale
  switch (scaleType) {

  case SteepestDescent::TwoNorm:

    dir.scale(-1.0/dir.norm());
    break;

  case SteepestDescent::FunctionTwoNorm:

    dir.scale(-1.0/soln.getNormF());
    break;

  case SteepestDescent::QuadMin:
  {
    // If necessary, allocate space for tmpVecPtr
    if (tmpVecPtr == NULL) 
      tmpVecPtr = soln.getX().clone(NOX::ShapeCopy);

    // Create a local reference
    Abstract::Vector& tmpVec(*tmpVecPtr);
    
    // Compute denominator
    ok = soln.applyJacobian(dir, tmpVec);
    if (!ok) {
      cerr << "NOX::Direction::SteepestDescent::compute - Unable to apply Jacobian" << endl;
      return false;
    }

    // Scale
    dir.scale( -1.0 * dir.dot(dir) / tmpVec.dot(tmpVec) );

    break;
  }
  case SteepestDescent::None:

    dir.scale( -1.0 );
    break;

  default:

    cout << "NOX::Direction::compute - Invalid scaleType" << endl;
    throw "NOX Error";

  }

  return true;
}


