// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
#include "LOCA_Continuation_NaturalGroup.H"
#include "LOCA_ErrorCheck.H"

LOCA::Continuation::NaturalGroup::NaturalGroup(
				  LOCA::Continuation::AbstractGroup& g,
				  int paramID,
				  NOX::Parameter::List& params)
  : LOCA::Continuation::ExtendedGroup(g, paramID, params),
    xVec(g.getX(), g.getParam(paramID)),
    fVec(g.getX(), 0.0),
    newtonVec(g.getX(), 0.0),
    gradientVec(g.getX(), 0.0),
    prevXVec(g.getX(), g.getParam(paramID)),
    derivResidualParamPtr(g.getX().clone(NOX::ShapeCopy)),
    stepSize(0.0),
    isValidPrevXVec(false)
{
  resetIsValid();
}

LOCA::Continuation::NaturalGroup::NaturalGroup(
				 LOCA::Continuation::AbstractGroup& g,
				 string paramID,
				 NOX::Parameter::List& params)
  : LOCA::Continuation::ExtendedGroup(g, paramID, params),
    xVec(g.getX(), g.getParam(paramID)),
    fVec(g.getX(), 0.0),
    newtonVec(g.getX(), 0.0),
    gradientVec(g.getX(), 0.0),
    prevXVec(g.getX(), g.getParam(paramID)),
    derivResidualParamPtr(g.getX().clone(NOX::ShapeCopy)),
    stepSize(0.0),
    isValidPrevXVec(false)
{
  resetIsValid();
}

LOCA::Continuation::NaturalGroup::NaturalGroup(
			     const LOCA::Continuation::NaturalGroup& source, 
			     NOX::CopyType type)
  : LOCA::Continuation::ExtendedGroup(source, type), 
    xVec(source.xVec, type),
    fVec(source.fVec, type),
    newtonVec(source.newtonVec, type),
    gradientVec(source.gradientVec, type),
    prevXVec(source.prevXVec, type),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    stepSize(source.stepSize),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    isValidPrevXVec(source.isValidPrevXVec)
{
}


LOCA::Continuation::NaturalGroup::~NaturalGroup() 
{
  delete derivResidualParamPtr;
}

LOCA::Continuation::ExtendedGroup&
LOCA::Continuation::NaturalGroup::operator=(
			     const LOCA::Continuation::ExtendedGroup& source)
{
  return *this = dynamic_cast<const LOCA::Continuation::NaturalGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Continuation::NaturalGroup::operator=(
				    const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const LOCA::Continuation::NaturalGroup&>(source);
}

LOCA::Continuation::NaturalGroup&
LOCA::Continuation::NaturalGroup::operator=(
			     const LOCA::Continuation::NaturalGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    LOCA::Continuation::ExtendedGroup::operator=(source);

    // Copy values
    xVec = source.xVec;
    fVec = source.fVec;
    newtonVec = source.newtonVec;
    gradientVec = source.gradientVec;
    prevXVec = source.prevXVec;
    *derivResidualParamPtr = *source.derivResidualParamPtr;
    stepSize = source.stepSize;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidPrevXVec = source.isValidPrevXVec;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Continuation::NaturalGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Continuation::NaturalGroup(*this, type);
}

void
LOCA::Continuation::NaturalGroup::setStepSize(double deltaS) 
{
  stepSize = deltaS;
  
  // stepsize appears on RHS but not in Jacobian
  isValidF = false;
  isValidNewton = false;
  isValidGradient = false;
}

double
LOCA::Continuation::NaturalGroup::getStepSize() const 
{
  return stepSize;
}

void 
LOCA::Continuation::NaturalGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y) );
}

void
LOCA::Continuation::NaturalGroup::setX(
				const LOCA::Continuation::ExtendedVector& y) 
{
  grpPtr->setX( y.getXVec() );
  grpPtr->setParam(conParamID, y.getParam());
  xVec = y;

  resetIsValid();
}

void
LOCA::Continuation::NaturalGroup::setPrevX(const NOX::Abstract::Vector& y) 
{
  setPrevX( dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y) );
}

void
LOCA::Continuation::NaturalGroup::setPrevX(
				 const LOCA::Continuation::ExtendedVector& y) 
{
  prevXVec = y;

  isValidF = false; // Previous vector part of residual equation
  isValidNewton = false;
  isValidGradient = false;
  isValidPrevXVec = true;
}

const LOCA::Continuation::ExtendedVector&
LOCA::Continuation::NaturalGroup::getPrevX() const 
{
  return prevXVec;
}

bool
LOCA::Continuation::NaturalGroup::isPrevXVec() const
{
  return isValidPrevXVec;
}

void
LOCA::Continuation::NaturalGroup::scalePredictor(
					LOCA::Continuation::ExtendedVector& v)
{
}

void
LOCA::Continuation::NaturalGroup::computeX(const NOX::Abstract::Group& g, 
					   const NOX::Abstract::Vector& d,
					   double step) 
{
  computeX( dynamic_cast<const LOCA::Continuation::NaturalGroup&>(g),
	    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(d),
	    step);
}

void
LOCA::Continuation::NaturalGroup::computeX(
				 const LOCA::Continuation::NaturalGroup& g, 
				 const LOCA::Continuation::ExtendedVector& d,
				 double step) 
{
  grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  xVec.update(1.0, g.getX(), step, d, 0.0);
  grpPtr->setParam(conParamID, xVec.getParam());

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  string callingFunction = "LOCA::Continuation::NaturalGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute underlying F
  if (!grpPtr->isF()) {
    finalStatus = grpPtr->computeF();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  fVec.getXVec() = grpPtr->getF();

  fVec.getParam() = xVec.getParam() - prevXVec.getParam() - stepSize;
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Continuation::NaturalGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
//   status = grpPtr->computeDfDp(conParamID, *derivResidualParamPtr);
//   finalStatus = 
//     LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						 callingFunction);

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Continuation::NaturalGroup::computeGradient()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  
  // Make sure Jacobian is valid
  if (!isJacobian()) {
    status = computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  
  // Compute underlying gradient
  if (!grpPtr->isGradient()) {
    status = grpPtr->computeGradient();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  gradientVec.getXVec() = grpPtr->getGradient();
  gradientVec.getParam() = 
    derivResidualParamPtr->innerProduct(fVec.getXVec())  + fVec.getParam();

  isValidGradient = true;

  return finalStatus;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeNewton(NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Continuation::NaturalGroup::computeNewton()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  
  // Make sure Jacobian is valid
  if (!isJacobian()) {
    status = computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // zero out newton vec -- used as initial guess for some linear solvers
  newtonVec.init(0.0);

  status = applyJacobianInverse(params, fVec, newtonVec);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  newtonVec.scale(-1.0);

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::applyJacobian(
				       const NOX::Abstract::Vector& input, 
				       NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Continuation::NaturalGroup::applyJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation vectors
  const LOCA::Continuation::ExtendedVector& c_input = 
    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(input);
  LOCA::Continuation::ExtendedVector& c_result = 
    dynamic_cast<LOCA::Continuation::ExtendedVector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::Vector& input_x = c_input.getXVec();
  double input_param = c_input.getParam();

  // Get references to x, param components of result vector
  NOX::Abstract::Vector& result_x = c_result.getXVec();
  double& result_param = c_result.getParam();

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
 
  // Parameter equation
  result_param = input_param;

  // compute J*x
  status = grpPtr->applyJacobian(input_x, result_x);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute J*x + p*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, 1.0);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::applyJacobianTranspose(
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Continuation::NaturalGroup::applyJacobianTranspose()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation vectors
  const LOCA::Continuation::ExtendedVector& c_input = 
    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(input);
  LOCA::Continuation::ExtendedVector& c_result = 
    dynamic_cast<LOCA::Continuation::ExtendedVector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::Vector& input_x = c_input.getXVec();
  double input_param = c_input.getParam();

  // Get references to x, param components of result vector
  NOX::Abstract::Vector& result_x = c_result.getXVec();
  double& result_param = c_result.getParam();
  
  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // compute J^T*x
  status = grpPtr->applyJacobianTranspose(input_x, result_x);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute df/dp^T x
  result_param = derivResidualParamPtr->innerProduct(input_x) + input_param;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::applyJacobianInverse(
					 NOX::Parameter::List& params, 
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Continuation::NaturalGroup::applyJacobianInverse()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation vectors
  const LOCA::Continuation::ExtendedVector& c_input = 
    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(input);
  LOCA::Continuation::ExtendedVector& c_result = 
    dynamic_cast<LOCA::Continuation::ExtendedVector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::Vector& input_x = c_input.getXVec();
  double input_param = c_input.getParam();

  // Get references to x, param components of result vector
  NOX::Abstract::Vector& result_x = c_result.getXVec();
  double& result_param = c_result.getParam();

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Parameter equation
  result_param = input_param;

  // Compute input_x + result_param*df/dp
  NOX::Abstract::Vector* a = input_x.clone(NOX::DeepCopy);
  //a->update(result_param, *derivResidualParamPtr, 1.0);

  // solve J*result_x = a
  status = grpPtr->applyJacobianInverse(params, *a, result_x);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  delete a;

  return finalStatus;
}

bool
LOCA::Continuation::NaturalGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Continuation::NaturalGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Continuation::NaturalGroup::isGradient() const 
{
  return isValidGradient;
}

bool
LOCA::Continuation::NaturalGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Continuation::NaturalGroup::getX() const 
{
  return xVec;
}

const NOX::Abstract::Vector&
LOCA::Continuation::NaturalGroup::getF() const 
{
  return fVec;
}

double
LOCA::Continuation::NaturalGroup::getNormF() const 
{
  return fVec.norm();
}

const NOX::Abstract::Vector&
LOCA::Continuation::NaturalGroup::getGradient() const 
{
  return gradientVec;
}

const NOX::Abstract::Vector&
LOCA::Continuation::NaturalGroup::getNewton() const 
{
  return newtonVec;
}

double
LOCA::Continuation::NaturalGroup::getNormNewtonSolveResidual() const 
{
  string callingFunction = 
    "LOCA::Continuation::NaturalGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::Continuation::ExtendedVector residual = fVec;
  
  finalStatus = applyJacobian(newtonVec, residual);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, fVec, 1.0);
  return residual.norm();
}

void
LOCA::Continuation::NaturalGroup::setContinuationParameter(double val) 
{
  LOCA::Continuation::ExtendedGroup::setContinuationParameter(val);
  xVec.getParam() = val;
  resetIsValid();
}

void
LOCA::Continuation::NaturalGroup::resetIsValid() {
  LOCA::Continuation::ExtendedGroup::resetIsValid();
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
}

