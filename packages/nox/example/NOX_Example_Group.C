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
#include "NOX_Example_Group.H"	// class definition

extern "C"
{
  void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* idb, int* info);
};

using namespace NOX;
using namespace NOX::Example;

Group::Group(Interface& interface):
  xVector(interface.getInitialGuess()),	// deep copy      
  fVector(xVector, ShapeCopy),	// new vector of same size
  gradientVector(xVector, ShapeCopy),   // new vector of same size
  newtonVector(xVector, ShapeCopy),	// new vector of same size
  jacobianMatrix(xVector.length()),	// create a Jacobian matrix
  problemInterface(interface)	        // set reference to the problem interface
{
  normF = 0;
  resetIsValid();
}

Group::Group(const Group& source, CopyType type) :
  xVector(source.xVector, type), 
  fVector(source.fVector, type), 
  gradientVector(source.gradientVector, type), 
  newtonVector(source.newtonVector, type),
  jacobianMatrix(source.jacobianMatrix, type),    
  problemInterface(source.problemInterface)	
{
 
  switch (type) {
    
  case DeepCopy:
    
    isValidF = source.isValidF;
    isValidGradient = source.isValidGradient;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    normF = source.normF;
    break;

  case ShapeCopy:
    resetIsValid();
    normF = 0.0;
    break;

  default:
    cerr << "NOX:Example::Group - invalid CopyType for copy constructor." << endl;
    throw "NOX Example Error";
  }

}

Group::~Group() 
{
}

void Group::resetIsValid() //private
{
  isValidF = false;
  isValidJacobian = false;
  isValidGradient = false;
  isValidNewton = false;
}

Abstract::Group* Group::clone(CopyType type) const 
{
  Group* newgrp = new Group(*this, type);
  return newgrp;
}

Abstract::Group& Group::operator=(const Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group& Group::operator=(const Group& source)
{
  // Copy the xVector
  xVector = source.xVector;

  // Update the isValidVectors
  isValidF = source.isValidF;
  isValidGradient = source.isValidGradient;
  isValidNewton = source.isValidNewton;
  isValidJacobian = source.isValidJacobian;

  // Only copy vectors that are valid
  if (isValidF) {
    fVector = source.fVector;
    normF = source.normF;
  }

  if (isValidGradient)
    gradientVector = source.gradientVector;

  if (isValidNewton)
    newtonVector = source.newtonVector;

  if (isValidJacobian)
    jacobianMatrix = source.jacobianMatrix;
    
  return *this;
}

bool Group::setX(const Abstract::Vector& y) 
{
  return setX(dynamic_cast<const Vector&> (y));
}

bool Group::setX(const Vector& y) 
{
  resetIsValid();
  xVector = y;
  return true;
}

bool Group::computeX(const Abstract::Group& grp, 
		     const Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& examplegrp = dynamic_cast<const Group&> (grp);
  const Vector& exampled = dynamic_cast<const Vector&> (d);
  return computeX(examplegrp, exampled, step); 
}

bool Group::computeX(const Group& grp, const Vector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return true;
}

bool Group::computeF() 
{
  if (isValidF)
    return true;

  isValidF = problemInterface.computeF(fVector, xVector);

  if (isValidF) 
    normF = fVector.norm();

  return isValidF;
}

bool Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isValidJacobian)
    return true;

  isValidJacobian = problemInterface.computeJacobian(jacobianMatrix, xVector);

  return isValidJacobian;
}

bool Group::computeGradient() 
{
  if (isValidGradient)
    return true;
  
  if (!isF()) {
    cerr << "ERROR: NOX::Example::Group::computeGrad() - F is out of date wrt X!" << endl;
    return false;
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Example::Group::computeGrad() - Jacobian is out of date wrt X!" << endl;
    return false;
  }
  
  // Compute Gradient = J' * F
  int n = xVector.length();
  for (int i = 0; i < n; i ++) { 
    gradientVector(i) = 0;
    for (int j = 0; j < n; j ++) {
      gradientVector(i) += jacobianMatrix(j,i) * xVector(j);
    }
  }

  // Reset isValidGradient
  isValidGradient = true;

  // Return result
  return true;
}

bool Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return true;

  if (!isF()) {
    cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid F" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid Jacobian" << endl;
    throw "NOX Error";
  }

  // Compute Newton Vector. Must copy jacobianMatrix and xVector into
  // temporary locations because the dgesv call overwrites them.
  int one = 1;
  int n = newtonVector.length();
  Matrix J(jacobianMatrix);
  Vector y(fVector);
  vector<int> ipiv(n,0);
  int info;

  ::dgesv_(&n, &one, &J(0,0), &n, &ipiv[0], &y(0), &n, &info);

  // Copy result into newtonVector
  newtonVector = y;

  // Scale soln by -1
  newtonVector.scale(-1.0);

  // Update state
  isValidNewton = (info == 0);

  // Return solution
  return isValidNewton;
}

bool Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& exampleinput = dynamic_cast<const Vector&> (input);
  Vector& exampleresult = dynamic_cast<Vector&> (result);
  return applyJacobian(exampleinput, exampleresult);
}

bool Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Compute result = J * input
  int n = input.length();
  for (int i = 0; i < n; i ++) { 
    result(i) = 0;
    for (int j = 0; j < n; j ++) {
      result(i) += jacobianMatrix(i,j) * input(j);
    }
  }

  return true;
}

bool Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& exampleinput = dynamic_cast<const Vector&> (input);
  Vector& exampleresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(exampleinput, exampleresult);
}

bool Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Compute result = J * input
  int n = input.length();
  for (int i = 0; i < n; i ++) { 
    result(i) = 0;
    for (int j = 0; j < n; j ++) {
      result(i) += jacobianMatrix(j,i) * input(j);
    }
  }

  return true;
}

bool Group::applyJacobianDiagonalInverse(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& exampleinput = dynamic_cast<const Vector&> (input);
  Vector& exampleresult = dynamic_cast<Vector&> (result);
  return applyJacobianDiagonalInverse(exampleinput, exampleresult);
}

bool Group::applyJacobianDiagonalInverse(const Vector& input, Vector& result) const
{
  if (!isJacobian()) 
    return false;

  // Compute result = J * input
  int n = input.length();
  for (int i = 0; i < n; i ++) { 
    result(i) = input(i) / jacobianMatrix(i,i);
  }

  return true;
}

bool Group::preconditionVector(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& exampleinput = dynamic_cast<const Vector&> (input);
  Vector& exampleresult = dynamic_cast<Vector&> (result);
  return preconditionVector(exampleinput, exampleresult);
}

bool Group::preconditionVector(const Vector& input, Vector& result) const
{
  cerr << "Warning - NOX::Example::Group::preconditionVector - not supported" << endl;
  return false;
}

bool Group::isF() const 
{   
  return isValidF;
}

bool Group::isJacobian() const 
{  
  return isValidJacobian;
}

bool Group::isGradient() const 
{   
  return isValidGradient;
}

bool Group::isNewton() const 
{   
  return isValidNewton;
}

const Abstract::Vector& Group::getX() const 
{
  return xVector;
}

const Abstract::Vector& Group::getF() const 
{  
  if (!isF()) {
    cerr << "ERROR: NOX::Example::Group::getF() - invalid F" << endl;
    throw "NOX Error";
  }
    
  return fVector;
}

double Group::getNormF() const
{
  return normF;
}

const Abstract::Vector& Group::getGradient() const 
{ 
  return gradientVector;
}

const Abstract::Vector& Group::getNewton() const 
{
  return newtonVector;
}


void Group::print() const
{
  cout << "x = " << xVector << "\n";

  if (isValidF) {
    cout << "F(x) = " << fVector << "\n";
    cout << "|| F(x) || = " << normF << "\n";
  }
  else
    cout << "F(x) has not been computed" << "\n";

  cout << endl;
}
