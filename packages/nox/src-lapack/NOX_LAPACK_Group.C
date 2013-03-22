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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Common.H"
#include "NOX_LAPACK_Group.H"	// class definition
#include "NOX_Abstract_MultiVector.H"
#include "Teuchos_BLAS_wrappers.hpp"
#include "Teuchos_LAPACK_wrappers.hpp"

NOX::LAPACK::Group::Group(NOX::LAPACK::Interface& interface):
  xVector(interface.getInitialGuess()),	// deep copy      
  fVector(xVector, ShapeCopy),	// new vector of same size
  newtonVector(xVector, ShapeCopy),	// new vector of same size
  gradientVector(xVector, ShapeCopy),   // new vector of same size
  jacSolver(xVector.length()),	// create a Jacobian matrix
  problemInterface(interface) // set reference to the problem interface
{
  resetIsValid();
}

NOX::LAPACK::Group::Group(const NOX::LAPACK::Group& source, NOX::CopyType type) :
  xVector(source.xVector, type), 
  fVector(source.fVector, type),  
  newtonVector(source.newtonVector, type),
  gradientVector(source.gradientVector, type),
  jacSolver(source.jacSolver),   
  problemInterface(source.problemInterface)	
{
 
  switch (type) {
    
  case NOX::DeepCopy:
    
    isValidF = source.isValidF;
    isValidGradient = source.isValidGradient;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    break;

  case NOX::ShapeCopy:
    resetIsValid();
    break;

  default:
    std::cerr << "NOX:LAPACK::Group - invalid CopyType for copy constructor." << std::endl;
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
  jacSolver.reset(); // Reset factorization
}

Teuchos::RCP<NOX::Abstract::Group> NOX::LAPACK::Group::
clone(NOX::CopyType type) const 
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp = 
    Teuchos::rcp(new NOX::LAPACK::Group(*this, type));
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
    }

    if (isValidGradient)
      gradientVector = source.gradientVector;
    
    if (isValidNewton)
      newtonVector = source.newtonVector;
    
    if (isValidJacobian)
      jacSolver = source.jacSolver;

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

  return (isValidF) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType NOX::LAPACK::Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  isValidJacobian = 
    problemInterface.computeJacobian(jacSolver.getMatrix(), xVector);

  return (isValidJacobian) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType NOX::LAPACK::Group::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;
  
  if (!isF()) {
    std::cerr << "ERROR: NOX::LAPACK::Group::computeGrad() - F is out of date wrt X!" << std::endl;
    return NOX::Abstract::Group::BadDependency;
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::LAPACK::Group::computeGrad() - Jacobian is out of date wrt X!" << std::endl;
    return NOX::Abstract::Group::BadDependency;
  }
  
  // Compute Gradient = J' * F
  jacSolver.apply(true, 1, &fVector(0), &gradientVector(0));

  // Reset isValidGradient
  isValidGradient = true;

  // Return result
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType NOX::LAPACK::Group::computeNewton(Teuchos::ParameterList& p) 
{
  if (isNewton())
    return NOX::Abstract::Group::Ok;

  if (!isF()) {
    std::cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid F" << std::endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid Jacobian" << std::endl;
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

  // Apply Jacobian
  jacSolver.apply(false, 1, &input(0), &result(0));

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

  // Apply Jacobian transpose
  jacSolver.apply(true, 1, &input(0), &result(0));

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const 
{
  const Vector& lapackinput = dynamic_cast<const Vector&> (input);
  Vector& lapackresult = dynamic_cast<Vector&> (result); 
  return applyJacobianInverse(p, lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const Vector& input, 
					 Vector& result) const 
{

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::LAPACK::Group::applyJacobianInverse() - invalid Jacobian" << std::endl;
    throw "NOX Error";
  }

  // Solve Jacobian
  result = input;
  bool res = jacSolver.solve(false, 1, &result(0));
    
  return res ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType 
NOX::LAPACK::Group::applyJacobianInverseMultiVector(
				     Teuchos::ParameterList& p, 
				     const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& result) const 
{

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::LAPACK::Group::applyJacobianInverseMultiVector() "
	 << "- invalid Jacobian" << std::endl;
    throw "NOX Error";
  }

  // Number of RHS
  int nVecs = input.numVectors();

  int m = jacSolver.getMatrix().numRows();

  // Copy all input vectors into one matrix
  NOX::LAPACK::Matrix<double> B(m,nVecs);
  const NOX::LAPACK::Vector* constVecPtr;
  for (int j=0; j<nVecs; j++) {
    constVecPtr = dynamic_cast<const NOX::LAPACK::Vector*>(&(input[j]));
    for (int i=0; i<m; i++)
      B(i,j) = (*constVecPtr)(i);
  }

  // Solve Jacobian
  bool res = jacSolver.solve(false, nVecs, &B(0,0));

  if (!res)
    return NOX::Abstract::Group::Failed;

  // Copy result from matrix
  NOX::LAPACK::Vector* vecPtr;
  for (int j=0; j<nVecs; j++) {
    vecPtr = dynamic_cast<NOX::LAPACK::Vector*>(&(result[j]));
    for (int i=0; i<m; i++)
      (*vecPtr)(i) = B(i,j);
  }
    
  return NOX::Abstract::Group::Ok;
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
  if (isValidF) 
    return fVector.norm();
  
  std::cerr << "ERROR: NOX::LAPACK::Group::getNormF() "
       << "- invalid F, please call computeF() first." << std::endl;
  throw "NOX Error";
  
  return 0.0;
}

const NOX::Abstract::Vector& NOX::LAPACK::Group::getGradient() const 
{ 
  return gradientVector;
}

const NOX::Abstract::Vector& NOX::LAPACK::Group::getNewton() const 
{
  return newtonVector;
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::LAPACK::Group::getXPtr() const 
{
  return Teuchos::rcp< const NOX::Abstract::Vector >(&xVector, false);
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::LAPACK::Group::getFPtr() const 
{  
  return Teuchos::rcp< const NOX::Abstract::Vector >(&fVector, false);
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::LAPACK::Group::getGradientPtr() const 
{ 
  return Teuchos::rcp< const NOX::Abstract::Vector >(&gradientVector, false);
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::LAPACK::Group::getNewtonPtr() const 
{
  return Teuchos::rcp< const NOX::Abstract::Vector >(&newtonVector, false);
}


void NOX::LAPACK::Group::print() const
{
  std::cout << "x = " << xVector << "\n";

  if (isValidF) {
    std::cout << "F(x) = " << fVector << "\n";
  }
  else
    std::cout << "F(x) has not been computed" << "\n";
  
  std::cout << std::endl;
}
