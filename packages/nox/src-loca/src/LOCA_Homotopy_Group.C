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

#include "LOCA_Homotopy_Group.H"      // class definition
#include "NOX_Parameter_List.H"      // data member
#include "LOCA_Parameter_Vector.H"

LOCA::Homotopy::Group::Group(NOX::Parameter::List& locaSublist,
			     LOCA::Homotopy::AbstractGroup& g,
			     double scalarRandom,
			     double scalarInitialGuess) :
  grpPtr(&g),
  gVecPtr(g.getX().clone(NOX::ShapeCopy)),
  randomVecPtr(gVecPtr->clone(NOX::ShapeCopy)),
  newtonVecPtr(0),
  gradVecPtr(0),
  paramVec(grpPtr->getParams()),
  conParam(0.0),
  conParamID(-1),
  conParamLabel("Homotopy Continuation Parameter"),
  ownsGroup(false)
{
  // construct a random vector for the problem 
  randomVecPtr->random();
  randomVecPtr->abs(*randomVecPtr);
  randomVecPtr->update(scalarInitialGuess, grpPtr->getX(), scalarRandom);

  // Set the isValid flags to false
  resetIsValidFlags();

  // Set the homotopy parameter as a parameter in the ParameterVector.
  // This will allow us to do an invasive homotopy since the app can now 
  // get access to the homotopy continuation parameter.
  paramVec.addParameter(conParamLabel, conParam);
  grpPtr->setParams(paramVec);

  // Set up paramID
  conParamID = paramVec.getIndex(conParamLabel);

  setStepperParameters(locaSublist);

}

LOCA::Homotopy::Group::Group(NOX::Parameter::List& locaSublist,
			     LOCA::Homotopy::AbstractGroup& g,
			     const NOX::Abstract::Vector& randomVector) :
  grpPtr(&g),
  gVecPtr(g.getX().clone(NOX::ShapeCopy)),
  randomVecPtr(gVecPtr->clone(NOX::ShapeCopy)),
  newtonVecPtr(0),
  gradVecPtr(0),
  paramVec(grpPtr->getParams()),
  conParam(0.0),
  conParamID(-1),
  conParamLabel("Homotopy Continuation Parameter"),
  ownsGroup(false)
{
  // construct a random vector for the problem 
  *randomVecPtr = randomVector;

  // Set the isValid flags to false
  resetIsValidFlags();

  // Set the homotopy parameter as a parameter in the ParameterVector.
  // This will allow us to do an invasive homotopy since the app can now 
  // get access to the homotopy continuation parameter.
  paramVec.addParameter(conParamLabel, conParam);
  grpPtr->setParams(paramVec);

  // Set up paramID
  conParamID = paramVec.getIndex(conParamLabel);

  setStepperParameters(locaSublist);

}

LOCA::Homotopy::Group::Group(const LOCA::Homotopy::Group& source, 
			     NOX::CopyType type) : 
  grpPtr(dynamic_cast<LOCA::Homotopy::AbstractGroup*>(source.grpPtr->clone(type))),
  gVecPtr(source.gVecPtr->clone(type)),
  randomVecPtr(source.randomVecPtr->clone(NOX::DeepCopy)), //deep copy to keep the set of equations consistent
  newtonVecPtr(0),
  gradVecPtr(0),
  paramVec(source.paramVec),
  conParam(source.conParam),
  conParamID(source.conParamID),
  conParamLabel(source.conParamLabel),
  ownsGroup(true)
{
  if (source.newtonVecPtr != 0)
    newtonVecPtr = source.newtonVecPtr->clone(type);

  if (source.gradVecPtr != 0)
    newtonVecPtr = source.gradVecPtr->clone(type);

  if (type == NOX::DeepCopy) {
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;
  }
  else if (type == NOX::ShapeCopy) {
    resetIsValidFlags();
  }
  else {
    errorCheck.throwError("LOCA::Homotopy::Group::Group(copy ctor)",
		     "CopyType is invalid!");
  }

}


LOCA::Homotopy::Group::~Group() 
{
  if (ownsGroup) 
    delete grpPtr;
  delete gVecPtr;
  delete randomVecPtr;
  delete newtonVecPtr;
  delete gradVecPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Homotopy::Group::operator=(
			     const LOCA::Continuation::AbstractGroup& source)
{
  return *this = dynamic_cast<const LOCA::Homotopy::Group&>(source);
}

LOCA::Extended::AbstractGroup&
LOCA::Homotopy::Group::operator=(
			     const LOCA::Extended::AbstractGroup& source)
{
  return *this = dynamic_cast<const LOCA::Homotopy::Group&>(source);
}

NOX::Abstract::Group&
LOCA::Homotopy::Group::operator=(const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const LOCA::Homotopy::Group&>(source);
}

LOCA::Homotopy::Group&
LOCA::Homotopy::Group::operator=(const LOCA::Homotopy::Group& source) 
{

  // Protect against A = A
  if (this != &source) {

    conParam = source.conParam;
    *grpPtr = *(source.grpPtr);
    *gVecPtr = *(source.gVecPtr);
    *randomVecPtr = *(source.randomVecPtr);
    if (newtonVecPtr != 0)
      *newtonVecPtr = *(source.newtonVecPtr);
    if (gradVecPtr != 0)
      *gradVecPtr = *(source.gradVecPtr);
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Homotopy::Group::clone(NOX::CopyType type) const 
{
  return new LOCA::Homotopy::Group(*this, type);
}

void
LOCA::Homotopy::Group::setParams(const LOCA::ParameterVector& p) 
{
  resetIsValidFlags();
  grpPtr->setParams(p);
  p.getValue(conParamLabel);
}

const LOCA::ParameterVector&
LOCA::Homotopy::Group::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::Homotopy::Group::setParam(int paramID, double val)
{
  resetIsValidFlags();
  grpPtr->setParam(paramID, val);
  if (paramID == conParamID)
    conParam = val;
}

double
LOCA::Homotopy::Group::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

void
LOCA::Homotopy::Group::setParam(string paramID, double val)
{
  resetIsValidFlags();
  grpPtr->setParam(paramID, val);
  if (paramID == conParamLabel)
    conParam = val;
}

double
LOCA::Homotopy::Group::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeDfDp(int paramID, 
				   NOX::Abstract::Vector& result)
{
  NOX::Abstract::Group::ReturnType res;

  // Compute f
  res = grpPtr->computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // g = conParam * f(x) + ((1.0 - conParam) * (x - randomVec))
  // dg/dp = f(x) - (x - randomVec) where p = conParam
  result = grpPtr->getF();
  result.update(-1.0, grpPtr->getX(), 1.0, *randomVecPtr, 1.0); 

  return res;
}

void
LOCA::Homotopy::Group::setX(const NOX::Abstract::Vector& y) 
{
  resetIsValidFlags();
  grpPtr->setX(y);
}

void
LOCA::Homotopy::Group::computeX(const NOX::Abstract::Group& g, 
			      const NOX::Abstract::Vector& d,
			      double step) 
{
  computeX( dynamic_cast<const LOCA::Homotopy::Group&>(g), d, step);
}

void
LOCA::Homotopy::Group::computeX(const LOCA::Homotopy::Group& g, 
			      const NOX::Abstract::Vector& d,
			      double step) 
{
  resetIsValidFlags();
  grpPtr->computeX(*(g.grpPtr), d, step);
}

NOX::Abstract::Group::ReturnType LOCA::Homotopy::Group::computeF()
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType status;

  status = grpPtr->computeF();
  errorCheck.checkReturnType(status, ErrorCheck::PrintWarning,
			"LOCA::Homotopy::Group::computeF",
			"grpPtr->computeF() failed!");
  if (status != NOX::Abstract::Group::Ok)
    return status;
  
  // g = conParam * f(x) + ((1.0 - conParam) * (x - randomVec))
  *gVecPtr = grpPtr->getX();
  gVecPtr->update(-1.0, *randomVecPtr, 1.0); 
  gVecPtr->scale( 1.0 - conParam);
  gVecPtr->update(conParam, grpPtr->getF(), 1.0);
 

  isValidF = true;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType status = grpPtr->computeJacobian();
  errorCheck.checkReturnType(status, ErrorCheck::PrintWarning,
			"LOCA::Homotopy::Group::computeJacobian",
			"grpPtr->computeJacobian() failed!");
  if (status != NOX::Abstract::Group::Ok)
    return status;

  // Augment the group's Jacobian for homotopy
  grpPtr->augmentJacobianForHomotopy(conParam);

  isValidJacobian = true;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType status = computeF();
  if (status != NOX::Abstract::Group::Ok)
    return status;
  
  status = computeJacobian();
  if (status != NOX::Abstract::Group::Ok)
    return status;

  status = applyJacobianTranspose(*gVecPtr, *gradVecPtr);

  return status;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeNewton(NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;
  
  if (newtonVecPtr == 0)
    newtonVecPtr = gVecPtr->clone(NOX::ShapeCopy);
  
  NOX::Abstract::Group::ReturnType status = computeF();
  if (status != NOX::Abstract::Group::Ok)
    return status;
  
  status = computeJacobian();
  if (status != NOX::Abstract::Group::Ok)
    return status;

  status = applyJacobianInverse(params, *gVecPtr, *newtonVecPtr);
  if (status != NOX::Abstract::Group::Ok)
    return status;

  newtonVecPtr->scale(-1.0);
  
  isValidNewton = true;
  
  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobian(const NOX::Abstract::Vector& input,
				   NOX::Abstract::Vector& result) const 
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;

  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobian(input, result);

  errorCheck.checkReturnType(status, ErrorCheck::PrintWarning,
			"LOCA::Homotopy::Group::applyJacobian",
			"grpPtr->applyJacobian() failed!");
  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianTranspose(
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;

  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianTranspose(input, result);
  
  errorCheck.checkReturnType(status, ErrorCheck::PrintWarning,
			"LOCA::Homotopy::Group::applyJacobianTranspose",
			"grpPtr->applyJacobianTranspose() failed!");
  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianInverse(NOX::Parameter::List& params,
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianInverse(params, input, result);

  errorCheck.checkReturnType(status, ErrorCheck::PrintWarning,
			"LOCA::Homotopy::Group::applyJacobianInverse",
			"grpPtr->applyJacobianInverse() failed!");
  return status;
}

bool
LOCA::Homotopy::Group::isF() const 
{
  return isValidF;
}

bool
LOCA::Homotopy::Group::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Homotopy::Group::isNewton() const 
{
  return isValidNewton;
}
  
bool
LOCA::Homotopy::Group::isGradient() const 
{
  return isValidGradient;
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getX() const 
{
  return grpPtr->getX();
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getF() const 
{
  return *gVecPtr;
}

double
LOCA::Homotopy::Group::getNormF() const 
{
  return gVecPtr->norm();
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getGradient() const 
{
  if (gradVecPtr == 0) {
    errorCheck.throwError("LOCA::Homotopy::Group::getGradient", 
		     "gradVecPtr is NULL!");
  }
  return *gradVecPtr;
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getNewton() const 
{
  if (newtonVecPtr == 0) {
    errorCheck.throwError("LOCA::Homotopy::Group::getNewton", 
		     "newtonVecPtr is NULL!");
  }
  return *newtonVecPtr;
}

void
LOCA::Homotopy::Group::printSolution(const double conParam) const 
{
  cout << "LOCA::Homotopy::Group::printSolution\n";
  cout << "Nothing to print at this time!";
  return;
}

void
LOCA::Homotopy::Group::resetIsValidFlags()
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
  return;
}

void
LOCA::Homotopy::Group::setStepperParameters(NOX::Parameter::List& params)
{
  
  // Create the stepper sublist
  NOX::Parameter::List& stepperList = params.sublist("Stepper");
  stepperList.setParameter("Continuation Method", "Natural");
  stepperList.setParameter("Continuation Parameter", conParamLabel);
  stepperList.setParameter("Initial Value", 0.0);
  stepperList.setParameter("Max Value", 1.0);
  stepperList.setParameter("Min Value", -1.0);
  stepperList.setParameter("Max Steps", 50);

  // Create predictor sublist
  NOX::Parameter::List& predictorList = params.sublist("Predictor");
  predictorList.setParameter("Method", "Constant");

  // Create step size sublist
  NOX::Parameter::List& stepSizeList = params.sublist("Step Size");
  //stepSizeList.setParameter("Method", "Constant");
  stepSizeList.setParameter("Method", "Adaptive");
  stepSizeList.setParameter("Initial Step Size", 0.1);
  stepSizeList.setParameter("Min Step Size", 1.0e-2);
  stepSizeList.setParameter("Max Step Size", 1.0);
  stepSizeList.setParameter("Aggressiveness", 0.5);
  return;
}

const LOCA::Continuation::AbstractGroup&
LOCA::Homotopy::Group::getUnderlyingGroup() const
{
  return *grpPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Homotopy::Group::getUnderlyingGroup()
{
  return *grpPtr;
}
