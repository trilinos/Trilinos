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
#include "LOCA_Continuation_NaturalGroup.H"

LOCA::Continuation::NaturalGroup::NaturalGroup(const LOCA::Abstract::Group& g,
					       int paramID,
					       NOX::Parameter::List& linSolverParams,
					       NOX::Parameter::List& params)
  : LOCA::Continuation::Group(g, paramID, linSolverParams, params),
    xVec(g.getX(), g.getParam(paramID)),
    fVec(g.getX(), 0.0),
    newtonVec(g.getX(), 0.0),
    prevXVec(g.getX(), g.getParam(paramID)),
    derivResidualParamPtr(g.getX().clone(NOX::ShapeCopy)),
    stepSize(0.0),
    isValidPrevXVec(false)
{
  resetIsValid();
}

LOCA::Continuation::NaturalGroup::NaturalGroup(const LOCA::Abstract::Group& g,
					       string paramID,
					       NOX::Parameter::List& linSolverParams,
					       NOX::Parameter::List& params)
  : LOCA::Continuation::Group(g, paramID, linSolverParams, params),
    xVec(g.getX(), g.getParam(paramID)),
    fVec(g.getX(), 0.0),
    newtonVec(g.getX(), 0.0),
    prevXVec(g.getX(), g.getParam(paramID)),
    derivResidualParamPtr(g.getX().clone(NOX::ShapeCopy)),
    stepSize(0.0),
    isValidPrevXVec(false)
{
  resetIsValid();
}

LOCA::Continuation::NaturalGroup::NaturalGroup(const LOCA::Continuation::NaturalGroup& source, NOX::CopyType type)
  : LOCA::Continuation::Group(source, type), 
    xVec(source.xVec, type),
    fVec(source.fVec, type),
    newtonVec(source.newtonVec, type),
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

LOCA::Continuation::Group&
LOCA::Continuation::NaturalGroup::operator=(const LOCA::Continuation::Group& source)
{
  return *this = dynamic_cast<const LOCA::Continuation::NaturalGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Continuation::NaturalGroup::operator=(const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const LOCA::Continuation::NaturalGroup&>(source);
}

LOCA::Continuation::NaturalGroup&
LOCA::Continuation::NaturalGroup::operator=(const LOCA::Continuation::NaturalGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    LOCA::Continuation::Group::operator=(source);

    // delete old values
    delete derivResidualParamPtr;

    // Copy values
    xVec = source.xVec;
    fVec = source.fVec;
    newtonVec = source.newtonVec;
    prevXVec = source.prevXVec;
    derivResidualParamPtr = source.derivResidualParamPtr->clone(NOX::DeepCopy);
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
}

double
LOCA::Continuation::NaturalGroup::getStepSize() 
{
  return stepSize;
}

void 
LOCA::Continuation::NaturalGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Continuation::Vector&>(y) );
}

void
LOCA::Continuation::NaturalGroup::setX(const LOCA::Continuation::Vector& y) 
{
  LOCA::Continuation::Group::grpPtr->setX( y.getXVec() );
  LOCA::Continuation::Group::grpPtr->setParam(LOCA::Continuation::Group::conParamID, y.getParam());
  xVec = y;

  resetIsValid();
}

void
LOCA::Continuation::NaturalGroup::setPrevX(const NOX::Abstract::Vector& y) 
{
  setPrevX( dynamic_cast<const LOCA::Continuation::Vector&>(y) );
}

void
LOCA::Continuation::NaturalGroup::setPrevX(const LOCA::Continuation::Vector& y) 
{
  prevXVec = y;

  isValidF = false; // Previous vector part of residual equation
  isValidNewton = false;
  isValidPrevXVec = true;
}

const LOCA::Continuation::Vector&
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
LOCA::Continuation::NaturalGroup::computeX(const NOX::Abstract::Group& g, 
					   const NOX::Abstract::Vector& d,
					   double step) 
{
  computeX( dynamic_cast<const LOCA::Continuation::NaturalGroup&>(g),
	    dynamic_cast<const LOCA::Continuation::Vector&>(d),
	    step);
}

void
LOCA::Continuation::NaturalGroup::computeX(const LOCA::Continuation::NaturalGroup& g, 
					   const LOCA::Continuation::Vector& d,
					   double step) 
{
  LOCA::Continuation::Group::grpPtr->computeX(*(g.LOCA::Continuation::Group::grpPtr), d.getXVec(), step);
  xVec.update(1.0, g.getX(), step, d, 0.0);
  LOCA::Continuation::Group::grpPtr->setParam(LOCA::Continuation::Group::conParamID, xVec.getParam());

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  res = LOCA::Continuation::Group::grpPtr->computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  fVec.getXVec() = LOCA::Continuation::Group::grpPtr->getF();

  fVec.getParam() = xVec.getParam() - prevXVec.getParam() - stepSize;
  
  isValidF = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  res =  LOCA::Continuation::Group::grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = LOCA::Continuation::Group::grpPtr->computeDfDp(LOCA::Continuation::Group::conParamID, *derivResidualParamPtr);

  if (res != NOX::Abstract::Group::Ok)
    return res;

  isValidJacobian = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::computeNewton(NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = applyJacobianInverse(params, fVec, newtonVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  newtonVec.scale(-1.0);

  isValidNewton = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
  // Cast inputs to continuation vectors
  const LOCA::Continuation::Vector& c_input = 
    dynamic_cast<const LOCA::Continuation::Vector&>(input);
  LOCA::Continuation::Vector& c_result = 
    dynamic_cast<LOCA::Continuation::Vector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::Vector& input_x = c_input.getXVec();
  double input_param = c_input.getParam();

  // Get references to x, param components of result vector
  NOX::Abstract::Vector& result_x = c_result.getXVec();
  double& result_param = c_result.getParam();

  NOX::Abstract::Group::ReturnType res;
 
  // Parameter equation
  result_param = input_param;

  // compute J*x
  res = LOCA::Continuation::Group::grpPtr->applyJacobian(input_x, result_x);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*x + p*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, 1.0);

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::applyJacobianInverse(NOX::Parameter::List& params, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
  // Cast inputs to continuation vectors
  const LOCA::Continuation::Vector& c_input = 
    dynamic_cast<const LOCA::Continuation::Vector&>(input);
  LOCA::Continuation::Vector& c_result = 
    dynamic_cast<LOCA::Continuation::Vector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::Vector& input_x = c_input.getXVec();
  double input_param = c_input.getParam();

  // Get references to x, param components of result vector
  NOX::Abstract::Vector& result_x = c_result.getXVec();
  double& result_param = c_result.getParam();
  
  NOX::Abstract::Group::ReturnType res;

  // Parameter equation
  result_param = input_param;

  // If result_param = 0, just solve J*result_x = input_x
  if (result_param = 0.0) {
    res = LOCA::Continuation::Group::grpPtr->applyJacobianInverse(params,
								  input_x,
								  result_x);
  }
  else {
    // Compute input_x + result_param*df/dp
    NOX::Abstract::Vector* a = input_x.clone(NOX::DeepCopy);
    a->update(result_param, *derivResidualParamPtr, 1.0);

    // solve J*result_x = a
    res = LOCA::Continuation::Group::grpPtr->applyJacobianInverse(params,
								  *a,
								  result_x);

    delete a;
  }

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::NaturalGroup::applyRightPreconditioning(NOX::Parameter::List& params, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
  // Cast inputs to continuation vectors
  const LOCA::Continuation::Vector& c_input = 
    dynamic_cast<const LOCA::Continuation::Vector&>(input);
  LOCA::Continuation::Vector& c_result = 
    dynamic_cast<LOCA::Continuation::Vector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::Vector& input_x = c_input.getXVec();
  double input_param = c_input.getParam();

  // Get references to x, param components of result vector
  NOX::Abstract::Vector& result_x = c_result.getXVec();
  double& result_param = c_result.getParam();
  
  NOX::Abstract::Group::ReturnType res;

  // Parameter equation
  result_param = input_param;

  // If result_param = 0, just solve J*result_x = input_x
  if (result_param = 0.0) {
    res = LOCA::Continuation::Group::grpPtr->applyRightPreconditioning(params,
								       input_x,
								       result_x);
  }
  else {
    // Compute input_x + result_param*df/dp
    NOX::Abstract::Vector* a = input_x.clone(NOX::DeepCopy);
    a->update(result_param, *derivResidualParamPtr, 1.0);

    // solve J*result_x = a
    res = LOCA::Continuation::Group::grpPtr->applyRightPreconditioning(params,
								       *a,
								       result_x);
    delete a;
  }

  return res;
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
  return false;
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
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Continuation::NaturalGroup::getNewton() const 
{
  return newtonVec;
}

double
LOCA::Continuation::NaturalGroup::getNormNewtonSolveResidual() const 
{
  LOCA::Continuation::Vector residual = fVec;
  
  NOX::Abstract::Group::ReturnType res = applyJacobian(newtonVec, residual);
  if (res != NOX::Abstract::Group::Ok) {
    cout << "ERROR: applyJacobian() in getNormNewtonSolveResidual "
	 << " returned not ok" << endl;
    throw "LOCA Error";
    return 0.0;
  }

  residual = residual.update(1.0, fVec, 1.0);
  return residual.norm();
}

void
LOCA::Continuation::NaturalGroup::setContinuationParameter(double val) 
{
  LOCA::Continuation::Group::setContinuationParameter(val);
  xVec.getParam() = val;
  resetIsValid();
}

void
LOCA::Continuation::NaturalGroup::resetIsValid() {
  LOCA::Continuation::Group::resetIsValid();
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

