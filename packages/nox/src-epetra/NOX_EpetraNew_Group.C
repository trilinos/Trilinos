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

#include "NOX_EpetraNew_Group.H"	// class definition

#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"
#include "NOX_EpetraNew_Interface_Required.H"
#include "Epetra_Vector.h"

using namespace NOX;
using namespace NOX::EpetraNew;

Group::Group(NOX::Parameter::List& printParams,
	     NOX::EpetraNew::Interface::Required& i, 
	     NOX::Epetra::Vector& x):
  utils(printParams),
  xVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(DeepCopy))),
  xVector(*xVectorPtr),    
  RHSVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  NewtonVector(*NewtonVectorPtr),
  tmpVectorPtr(0),
  normRHS(0.0),
  normNewtonSolveResidual(0),
  sharedLinearSystemPtr(0), 
  sharedLinearSystem(*sharedLinearSystemPtr),
  userInterface(i)
{
  // Set all isValid flags to false
  resetIsValid();
}

Group::Group(NOX::Parameter::List& printParams,
		   NOX::EpetraNew::Interface::Required& i, 
		   NOX::Epetra::Vector& x,
		   NOX::EpetraNew::LinearSystem& linSys):
  utils(printParams),
  xVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(DeepCopy))),
  xVector(*xVectorPtr),    
  RHSVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(x.clone(ShapeCopy))),
  NewtonVector(*NewtonVectorPtr),
  tmpVectorPtr(0),
  normRHS(0.0),
  normNewtonSolveResidual(0),
  sharedLinearSystemPtr(new NOX::SharedObject<NOX::EpetraNew::LinearSystem, NOX::EpetraNew::Group>(linSys)), 
  sharedLinearSystem(*sharedLinearSystemPtr),
  userInterface(i)
{
  // Set all isValid flags to false
  resetIsValid();
}

Group::Group(const Group& source, CopyType type) :
  utils(source.utils),
  xVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.xVector.clone(type))),
  xVector(*xVectorPtr),    
  RHSVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.RHSVector.clone(type))),
  RHSVector(*RHSVectorPtr), 
  gradVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.gradVector.clone(type))),
  gradVector(*gradVectorPtr), 
  NewtonVectorPtr(dynamic_cast<NOX::Epetra::Vector*>(source.NewtonVector.clone(type))),
  NewtonVector(*NewtonVectorPtr), 
  tmpVectorPtr(0),
  sharedLinearSystemPtr(0),
  sharedLinearSystem(source.sharedLinearSystem),
  userInterface(source.userInterface)
{
 
  switch (type) {
    
  case DeepCopy:
    
    isValidRHS = source.isValidRHS;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
    normRHS = source.normRHS;
    normNewtonSolveResidual = source.normNewtonSolveResidual;
    isValidPreconditioner = source.isValidPreconditioner;
    
    // New copy takes ownership of the shared Jacobian
    if (isValidJacobian)
      sharedLinearSystem.getObject(this);

    break;

  case ShapeCopy:
    resetIsValid();
    break;

  default:
    cerr << "ERROR: Invalid ConstructorType for group copy constructor." << endl;
    throw "NOX Error";
  }

}

Group::~Group() 
{
  delete tmpVectorPtr;
  delete sharedLinearSystemPtr;
  delete NewtonVectorPtr;
  delete gradVectorPtr;
  delete RHSVectorPtr;
  delete xVectorPtr;
}

void Group::resetIsValid() //private
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
  isValidNormNewtonSolveResidual = false;
  isValidPreconditioner = false;
  return;
}

Abstract::Group* Group::clone(CopyType type) const 
{
  NOX::Abstract::Group* newgrp = new NOX::EpetraNew::Group(*this, type);
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

  // Copy reference to sharedJacobian
  //sharedLinearSystemPtr = source.sharedLinearSystemPtr;

  // Update the isValidVectors
  isValidRHS = source.isValidRHS;
  isValidGrad = source.isValidGrad;
  isValidNewton = source.isValidNewton;
  isValidJacobian = source.isValidJacobian;
  isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
  isValidPreconditioner = source.isValidPreconditioner;

  // Only copy vectors that are valid
  if (isValidRHS) {
    RHSVector = source.RHSVector;
    normRHS = source.normRHS;
  }

  if (isValidGrad)
    gradVector = source.gradVector;

  if (isValidNewton)
    NewtonVector = source.NewtonVector;

  if (isValidNormNewtonSolveResidual)
    normNewtonSolveResidual = source.normNewtonSolveResidual;

  // If valid, this takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedLinearSystem.getObject(this);
    
  return *this;
}

void Group::setX(const Abstract::Vector& y) 
{
  setX(dynamic_cast<const NOX::Epetra::Vector&> (y));
  return;
}

void Group::setX(const NOX::Epetra::Vector& y) 
{
  if (isPreconditioner()) {
    sharedLinearSystem.getObject(this).destroyPreconditioner();
  }
  resetIsValid();
  xVector = y;
  return;
}

void Group::computeX(const NOX::Abstract::Group& grp, 
		     const NOX::Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const NOX::EpetraNew::Group& epetragrp = dynamic_cast<const Group&> (grp);
  const NOX::Epetra::Vector& epetrad = 
    dynamic_cast<const NOX::Epetra::Vector&> (d);
  computeX(epetragrp, epetrad, step); 
  return;
}

void Group::computeX(const Group& grp, 
		     const NOX::Epetra::Vector& d, 
		     double step) 
{
  if (isPreconditioner())
    sharedLinearSystem.getObject(this).destroyPreconditioner();
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return;
}

Abstract::Group::ReturnType Group::computeF() 
{
  if (isF())
    return Abstract::Group::Ok;

  bool status = false;
  
  status = userInterface.computeF(xVector.getEpetraVector(), RHSVector.getEpetraVector(), NOX::EpetraNew::Interface::Required::Residual);

  if (status == false) {
    cout << "ERROR: EpetraNew::Group::computeF() - fill failed!!!"
	 << endl;
    throw "NOX Error: Fill Failed";
  } 

  normRHS = RHSVector.norm();

  isValidRHS = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return Abstract::Group::Ok;

  // Fill the Jacobian 
  bool status = false;
  
  status = sharedLinearSystem.getObject(this).
    computeJacobian(xVector.getEpetraVector());

  if (status == false) {
    cout << "ERROR: NOX::EpetraNew::Group::computeJacobian() - fill failed!!!"
	 << endl;
    throw "NOX Error: Fill Failed";
  } 

  // Update status of Jacobian wrt solution vector
  isValidJacobian = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeGradient() 
{
  if (isGradient())
    return Abstract::Group::Ok;
  
  if (!isF()) {
    cerr << "ERROR: NOX::EpetraNew::Group::computeGradient() - RHS is out of date wrt X!" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::EpetraNew::Group::computeGradient() - Jacobian is out of date wrt X!" << endl;
    throw "NOX Error";
  }
  
  // Compute grad = Jacobian^T * RHS.
  sharedLinearSystem.getObject(this).applyJacobianTranspose(RHSVector,
							    gradVector);

  // Update state
  isValidGrad = true;

  // Return result
  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return Abstract::Group::Ok;

  if (!isF()) {
    cerr << "ERROR: NOX::EpetraNew::Group::computeNewton() - invalid RHS" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::EpetraNew::Group::computeNewton() - invalid Jacobian" << endl;
    throw "NOX Error";
  }

  Abstract::Group::ReturnType status;
  
  // Zero out the Newton Vector
  NewtonVector.init(0.0);

  // Create Epetra problem for the linear solve
  status = applyJacobianInverse(p, RHSVector, NewtonVector);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state EVEN IF LINEAR SOLVE FAILED
  // We still may want to use the vector even it it just missed it's 
  isValidNewton = true;

  // Compute the 2-norm of the linear solve residual ||Js+f||
  computeNormNewtonSolveResidual();

  // Return solution
  return status;
}

Abstract::Group::ReturnType Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetrainput = 
    dynamic_cast<const NOX::Epetra::Vector&> (input);
  NOX::Epetra::Vector& epetraresult = 
    dynamic_cast<NOX::Epetra::Vector&> (result);
  return applyJacobian(epetrainput, epetraresult);
}

Abstract::Group::ReturnType Group::applyJacobian(const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;

  // Apply the Jacobian
  bool status = sharedLinearSystem.getObject().applyJacobian(input, result);

  return status == true ? Abstract::Group::Ok : Abstract::Group::Failed;
}

Abstract::Group::ReturnType Group::applyJacobianInverse (Parameter::List &p, const Abstract::Vector &input, Abstract::Vector &result) const
{
  const NOX::Epetra::Vector& epetraInput = dynamic_cast<const NOX::Epetra::Vector&>(input);
  NOX::Epetra::Vector& epetraResult = dynamic_cast<NOX::Epetra::Vector&>(result);
  return applyJacobianInverse(p, epetraInput, epetraResult);
}

Abstract::Group::ReturnType Group::applyJacobianInverse (Parameter::List &p, const NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const
{
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;

  bool reusePrec = 
    sharedLinearSystem.getObject(this).checkPreconditionerReuse();

  if (!isPreconditioner()  && !reusePrec ) {
    sharedLinearSystem.getObject(this).destroyPreconditioner();
    sharedLinearSystem.getObject(this)
      .createPreconditioner(xVector.getEpetraVector(), p, false);
    isValidPreconditioner = true;
  }

  bool status = sharedLinearSystem.getObject(this).applyJacobianInverse(p, input, result);

  return status == true ? Abstract::Group::Ok : Abstract::Group::NotConverged;
}

Abstract::Group::ReturnType Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetrainput = dynamic_cast<const NOX::Epetra::Vector&> (input);
  NOX::Epetra::Vector& epetraresult = dynamic_cast<NOX::Epetra::Vector&> (result);
  return applyJacobianTranspose(epetrainput, epetraresult);
}

Abstract::Group::ReturnType Group::applyJacobianTranspose(const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;
  
  bool status = sharedLinearSystem.getObject().applyJacobianTranspose(input, result);

  return status == true ? Abstract::Group::Ok : Abstract::Group::Failed;
}

Abstract::Group::ReturnType Group::applyRightPreconditioning(
				      bool useTranspose,
				      Parameter::List& params,
				      const Abstract::Vector& input, 
				      Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetraInput = dynamic_cast<const NOX::Epetra::Vector&>(input);
  NOX::Epetra::Vector& epetraResult = dynamic_cast<NOX::Epetra::Vector&>(result);
  
  return applyRightPreconditioning(useTranspose, params, epetraInput, epetraResult);
}

Abstract::Group::ReturnType Group::applyRightPreconditioning(
				       bool useTranspose, 
				       Parameter::List& linearSolverParams,
				       const NOX::Epetra::Vector& input, 
				       NOX::Epetra::Vector& result) const
{

  bool success = false;
  
  if (!isPreconditioner()) {
    sharedLinearSystem.getObject(this).destroyPreconditioner();
    sharedLinearSystem.getObject(this).createPreconditioner(xVector.getEpetraVector(), linearSolverParams, false);
    isValidPreconditioner = true;
  }

  success = sharedLinearSystem.getObject()
    .applyRightPreconditioning(useTranspose, linearSolverParams, input, result);

  if (success == true) 
    return Abstract::Group::Ok;
  else
    return Abstract::Group::Failed;
}

bool Group::isF() const 
{   
  return isValidRHS;
}

bool Group::isJacobian() const 
{  
  return ((sharedLinearSystem.isOwner(this)) && (isValidJacobian));
}

bool Group::isGradient() const 
{   
  return isValidGrad;
}

bool Group::isNewton() const 
{   
  return isValidNewton;
}

bool Group::isNormNewtonSolveResidual() const 
{   
  return isValidNormNewtonSolveResidual;
}

bool Group::isPreconditioner() const 
{  
  return ((sharedLinearSystem.isOwner(this)) && (isValidPreconditioner));
}

const Abstract::Vector& Group::getX() const 
{
  return xVector;
}

const Abstract::Vector& Group::getF() const 
{  
  if (!isF()) {
    cerr << "ERROR: NOX::EpetraNew::Group::getF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return RHSVector;
}

double Group::getNormF() const
{
  if (!isF()) {
    cerr << "ERROR: NOX::EpetraNew::Group::getNormF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return normRHS;
}

const Abstract::Vector& Group::getGradient() const 
{ 
  if (!isGradient()) {
    cerr << "ERROR: NOX::EpetraNew::Group::getGradient() - invalid gradient" << endl;
    throw "NOX Error";
  }
    
  return gradVector;
}

const Abstract::Vector& Group::getNewton() const 
{
  if (!isNewton()) {
    cerr << "ERROR: NOX::EpetraNew::Group::getNewton() - invalid Newton vector" << endl;
    throw "NOX Error";
  }
    
  return NewtonVector;
}

Abstract::Group::ReturnType NOX::EpetraNew::Group::getNormLastLinearSolveResidual(double& residual) const
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual) {
    residual = normNewtonSolveResidual;
    return NOX::Abstract::Group::Ok;
  }
  
  // Otherwise give warning since a Newton direction has not been calculated
  // wrt this solution group
  if (utils.isPrintProcessAndType(Utils::Warning)) {
    cout << "ERROR: NOX::EpetraNew::Group::getNormLastLinearSolveResidual() - "
	 << "Group has not performed a Newton solve corresponding to this "
	 << "solution vector!" << endl;
  }
  return NOX::Abstract::Group::BadDependency;
}  
  
NOX::EpetraNew::Interface::Required& Group::getRequiredInterface()
{
  return userInterface;
}
 
const NOX::EpetraNew::LinearSystem& Group::getLinearSystem() const
{
  return sharedLinearSystem.getObject();
}

bool Group::computeNormNewtonSolveResidual ()
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual) 
    return true;

  // Make sure NewtonVector and RHSVector are valid
  // We could return false, but for now we will throw errors
  if (!isValidRHS) {
    cerr << "ERROR: NOX::EpetraNew::Group::computeNormNewtonSolveResidual() - invalid RHS" 
	 << endl;
    throw "NOX Error";
  }
  if (!isValidNewton) {
    cerr << "ERROR: NOX::EpetraNew::Group::computeNormNewtonSolveResidual() - invalid "
	 << "Newton direction" << endl;
    throw "NOX Error";
  }
  
  // Allocate the tmpVectorPtr if not already done (deleted in ~Group)   
  if (tmpVectorPtr == 0) {
    tmpVectorPtr = new Epetra_Vector(RHSVector.getEpetraVector());
  }
  NOX::Epetra::Vector tmpNoxVector(*tmpVectorPtr, ShapeCopy); 

  sharedLinearSystem.getObject().applyJacobian(NewtonVector, tmpNoxVector);    
  tmpNoxVector.update(1.0, RHSVector, 1.0);
  normNewtonSolveResidual = tmpNoxVector.norm();
  
  isValidNormNewtonSolveResidual = true;
  
  return true;
}
