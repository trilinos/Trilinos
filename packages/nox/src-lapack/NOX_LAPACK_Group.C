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
#include "NOX_LAPACK_Group.H"	// class definition
#include "NOX_BLAS_Wrappers.H"
#include "NOX_LAPACK_Wrappers.H"

NOX::LAPACK::Group::Group(NOX::LAPACK::Interface& interface):
  xVector(interface.getInitialGuess()),	// deep copy      
  fVector(xVector, ShapeCopy),	// new vector of same size
  gradientVector(xVector, ShapeCopy),   // new vector of same size
  newtonVector(xVector, ShapeCopy),	// new vector of same size
  jacobianMatrix(xVector.length(), xVector.length()),	// create a Jacobian matrix
  problemInterface(interface) // set reference to the problem interface
{
  normF = 0;
  resetIsValid();
}

NOX::LAPACK::Group::Group(const NOX::LAPACK::Group& source, NOX::CopyType type) :
  xVector(source.xVector, type), 
  fVector(source.fVector, type), 
  gradientVector(source.gradientVector, type), 
  newtonVector(source.newtonVector, type),
  jacobianMatrix(source.jacobianMatrix, type),   
  problemInterface(source.problemInterface)	
{
 
  switch (type) {
    
  case NOX::DeepCopy:
    
    isValidF = source.isValidF;
    isValidGradient = source.isValidGradient;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    normF = source.normF;
    break;

  case NOX::ShapeCopy:
    resetIsValid();
    normF = 0.0;
    break;

  default:
    cerr << "NOX:LAPACK::Group - invalid CopyType for copy constructor." << endl;
    throw "NOX LAPACK Error";
  }

}

NOX::LAPACK::Group::~Group() 
{
}

void NOX::LAPACK::Group::resetIsValid() //private
{
  isValidF = false;
  isValidJacobian = false;
  isValidGradient = false;
  isValidNewton = false;
}

NOX::Abstract::Group* NOX::LAPACK::Group::clone(NOX::CopyType type) const 
{
  Group* newgrp = new Group(*this, type);
  return newgrp;
}

NOX::Abstract::Group& NOX::LAPACK::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

NOX::Abstract::Group& NOX::LAPACK::Group::operator=(const Group& source)
{
  if (this != &source) {

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

  }

  return *this;
}

void NOX::LAPACK::Group::setX(const NOX::Abstract::Vector& y) 
{
  setX(dynamic_cast<const Vector&> (y));
}

void NOX::LAPACK::Group::setX(const Vector& y) 
{
  resetIsValid();
  xVector = y;
}

void NOX::LAPACK::Group::computeX(const NOX::Abstract::Group& grp, 
		     const NOX::Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& lapackgrp = dynamic_cast<const Group&> (grp);
  const Vector& lapackd = dynamic_cast<const Vector&> (d);
  computeX(lapackgrp, lapackd, step); 
}

void NOX::LAPACK::Group::computeX(const Group& grp, const Vector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
}

NOX::Abstract::Group::ReturnType NOX::LAPACK::Group::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  isValidF = problemInterface.computeF(fVector, xVector);

  if (isValidF) 
    normF = fVector.norm();

  return (isValidF) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType NOX::LAPACK::Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  isValidJacobian = problemInterface.computeJacobian(jacobianMatrix, xVector);

  return (isValidJacobian) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType NOX::LAPACK::Group::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;
  
  if (!isF()) {
    cerr << "ERROR: NOX::LAPACK::Group::computeGrad() - F is out of date wrt X!" << endl;
    return NOX::Abstract::Group::BadDependency;
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::LAPACK::Group::computeGrad() - Jacobian is out of date wrt X!" << endl;
    return NOX::Abstract::Group::BadDependency;
  }
  
  // Compute Gradient = J' * F
  int n = xVector.length();
  DGEMV_F77("T", &n, &n, &d_one, &jacobianMatrix(0,0), &n, &fVector(0),
	    &i_one, &d_zero, &gradientVector(0), &i_one);

  // Reset isValidGradient
  isValidGradient = true;

  // Return result
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType NOX::LAPACK::Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return NOX::Abstract::Group::Ok;

  if (!isF()) {
    cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid F" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid Jacobian" << endl;
    throw "NOX Error";
  }

  NOX::Abstract::Group::ReturnType status = applyJacobianInverse(p, fVector, newtonVector);
  isValidNewton = (status == NOX::Abstract::Group::Ok);

  // Scale soln by -1
  newtonVector.scale(-1.0);

  // Return solution
  return status;
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobian(const Abstract::Vector& input, 
				  NOX::Abstract::Vector& result) const
{
  const Vector& lapackinput = dynamic_cast<const Vector&> (input);
  Vector& lapackresult = dynamic_cast<Vector&> (result);
  return applyJacobian(lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return NOX::Abstract::Group::BadDependency;

  // Compute result = J * input
  int n = input.length();

  DGEMV_F77("N", &n, &n, &d_one, &jacobianMatrix(0,0), &n, &input(0),
	    &i_one, &d_zero, &result(0), &i_one);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobianTranspose(const NOX::Abstract::Vector& input, 
					   NOX::Abstract::Vector& result) const
{
  const Vector& lapackinput = dynamic_cast<const Vector&> (input);
  Vector& lapackresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return NOX::Abstract::Group::BadDependency;

  // Compute result = J * input
  int n = input.length();
  DGEMV_F77("T", &n, &n, &d_one, &jacobianMatrix(0,0), &n, &input(0),
	    &i_one, &d_zero, &result(0), &i_one);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobianInverse(NOX::Parameter::List& p, 
					 const Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const 
{
  const Vector& lapackinput = dynamic_cast<const Vector&> (input);
  Vector& lapackresult = dynamic_cast<Vector&> (result); 
  return applyJacobianInverse(p, lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobianInverse(NOX::Parameter::List& p, 
					 const Vector& input, Vector& result) const 
{

  if (!isF()) {
    cerr << "ERROR: NOX::LAPACK::Group::applyJacobianInverse() - invalid F" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::LAPACK::Group::applyJacobianInverse() - invalid Jacobian" << endl;
    throw "NOX Error";
  }

  // Must copy jacobianMatrix and xVector into temporary locations 
  // because the dgesv call overwrites them.
  int n = input.length();
  int info;
  Matrix J(jacobianMatrix);
  vector<int> ipiv(n,0);

  result = input;
  DGESV_F77(&n, &i_one, &J(0,0), &n, &ipiv[0], &result(0), &n, &info);
    
  return (info == 0) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

bool NOX::LAPACK::Group::isF() const 
{   
  return isValidF;
}

bool NOX::LAPACK::Group::isJacobian() const 
{  
  return isValidJacobian;
}

bool NOX::LAPACK::Group::isGradient() const 
{   
  return isValidGradient;
}

bool NOX::LAPACK::Group::isNewton() const 
{   
  return isValidNewton;
}

const NOX::Abstract::Vector& NOX::LAPACK::Group::getX() const 
{
  return xVector;
}

const NOX::Abstract::Vector& NOX::LAPACK::Group::getF() const 
{  
  return fVector;
}

double NOX::LAPACK::Group::getNormF() const
{
  return normF;
}

const NOX::Abstract::Vector& NOX::LAPACK::Group::getGradient() const 
{ 
  return gradientVector;
}

const NOX::Abstract::Vector& NOX::LAPACK::Group::getNewton() const 
{
  return newtonVector;
}


void NOX::LAPACK::Group::print() const
{
  cout << "x = " << xVector << "\n";

  if (isValidF) {
    cout << "F(x) = " << fVector << "\n";
    cout << "|| F(x) || = " << normF << "\n";
  }
  else
    cout << "F(x) has not been computed" << "\n";

  // if (isValidJacobian) {
  //     cout << "J(x) = " << endl << jacobianMatrix << "\n";
  //   }
  //   else
  //     cout << "J(x) has not been computed" << "\n";
  
  cout << endl;
}
