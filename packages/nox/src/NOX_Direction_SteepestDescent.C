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
#include "NOX_Utils.H"

NOX::Direction::SteepestDescent::SteepestDescent(NOX::Parameter::List& params) :
  tmpVecPtr(NULL)
{
  reset(params);
}

NOX::Direction::SteepestDescent::~SteepestDescent()
{
  delete tmpVecPtr;
}

bool NOX::Direction::SteepestDescent::reset(NOX::Parameter::List& params)
{
  const string tmp = params.getParameter("Scaling Type", "2-Norm");
  if (tmp == "2-Norm")
    scaleType = NOX::Direction::SteepestDescent::TwoNorm;
  else if (tmp == "F 2-Norm")
    scaleType = NOX::Direction::SteepestDescent::FunctionTwoNorm;
  else if (tmp == "Quadratic Model Min")
    scaleType = NOX::Direction::SteepestDescent::QuadMin;
  else if (tmp == "None")
    scaleType = NOX::Direction::SteepestDescent::None;
  else {
    cout << "NOX::Direction::SteepestDescent::reset - Invalid choice \""
	  << tmp << "\" for \"Scaling Type\"" << endl;
    throw "NOX Error";
  }

 return true;
}

bool NOX::Direction::SteepestDescent::compute(Abstract::Vector& dir, 
				 Abstract::Group& soln, 
				 const Solver::Generic& solver) 
{
  NOX::Abstract::Group::ReturnType status;


  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute F");

  // Compute Jacobian at current solution.
  status = soln.computeJacobian();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute Jacobian");

  // Compute the gradient at the current solution
  status = soln.computeGradient();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute gradient");

  // Get the gradient direction.
  dir = soln.getGradient();

  // Scale
  switch (scaleType) {

  case NOX::Direction::SteepestDescent::TwoNorm:

    dir.scale(-1.0/dir.norm());
    break;

  case NOX::Direction::SteepestDescent::FunctionTwoNorm:

    dir.scale(-1.0/soln.getNormF());
    break;

  case NOX::Direction::SteepestDescent::QuadMin:
  {
    // If necessary, allocate space for tmpVecPtr
    if (tmpVecPtr == NULL) 
      tmpVecPtr = soln.getX().clone(NOX::ShapeCopy);

    // Create a local reference
    Abstract::Vector& tmpVec(*tmpVecPtr);
    
    // Compute denominator
    status = soln.applyJacobian(dir, tmpVec);
    if (status != NOX::Abstract::Group::Ok) 
      throwError("compute", "Unable to compute apply Jacobian");

    // Scale
    dir.scale( -1.0 * dir.dot(dir) / tmpVec.dot(tmpVec) );

    break;
  }

  case NOX::Direction::SteepestDescent::None:

    dir.scale( -1.0 );
    break;

  default:
    
    throwError("compute", "Invalid scaleType");

  }

  return true;
}


void NOX::Direction::SteepestDescent::throwError(const string& functionName, const string& errorMsg)
{
    if (Utils::doPrint(Utils::Error))
      cerr << "NOX::Direction::SteepestDescent::" << functionName << " - " << errorMsg << endl;
    throw "NOX Error";
}
