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

#include "NOX_LineSearch_Utils_Slope.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"

NOX::LineSearch::Utils::Slope::Slope() :
  vecPtr(0)
{

}

NOX::LineSearch::Utils::Slope::~Slope()
{
  delete vecPtr;
}

double NOX::LineSearch::Utils::Slope::computeSlope(const Abstract::Vector& dir, const Abstract::Group& grp) 
{
   if (grp.isGradient()) 
     return(dir.dot(grp.getGradient()));

  // Allocate space for vecPtr if necessary
  if (vecPtr == 0) 
    vecPtr = dir.clone(ShapeCopy);

  // v = J * dir
  NOX::Abstract::Group::ReturnType status = grp.applyJacobian(dir,*vecPtr);
  
  if (status != NOX::Abstract::Group::Ok) 
  {
    cout << "NOX::LineSearch::Common::computeSlope -  Unable to apply Jacobian!" << endl;
    throw "NOX Error";
  }

  // Check that F exists
  if (!grp.isF()) 
  {
    cout << "NOX::LineSearch::Common::computeSlope - Invalid F" << endl;
    throw "NOX Error";
  }

  // Return <v, F> = F' * J * dir = <J'F, dir> = <g, dir>
  return(vecPtr->dot(grp.getF()));
}

double NOX::LineSearch::Utils::Slope::computeSlopeWithOutJac(const Abstract::Vector& dir, const Abstract::Group& grp) 
{
  // Allocate space for vecPtr and grpPtr if necessary
  if (vecPtr == 0) 
    vecPtr = dir.clone(ShapeCopy);
  if (grpPtr == 0)
    grpPtr = grp.clone(ShapeCopy);

  // Check that F exists
  if (!grp.isF()) 
  {
    cout << "NOX::LineSearch::Common::computeSlope - Invalid F" << endl;
    throw "NOX Error";
  }

  // Compute the perturbation parameter
  double lambda = 1.0e-6;
  double denominator = dir.norm();

  // Don't divide by zero
  if (denominator == 0.0)
    denominator = 1.0;

  double eta = lambda * (lambda + grp.getX().norm() / denominator);

  // Don't divide by zero
  if (eta == 0.0)
    eta = 1.0e-6;

  // Perturb the solution vector
  *vecPtr = grp.getX();
  vecPtr->update(eta, dir, 1.0);

  // Compute the new F --> F(x + eta * dir)
  grpPtr->setX(*vecPtr);  
  grpPtr->computeF();

  // Compute Js = (F(x + eta * dir) - F(x))/eta
  *vecPtr = grpPtr->getF();
  vecPtr->update(-1.0, grp.getF(), 1.0);
  vecPtr->scale(1.0/eta);
  
  return(vecPtr->dot(grp.getF()));
}
