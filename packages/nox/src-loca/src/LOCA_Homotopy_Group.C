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

#include "LOCA_Homotopy_Group.H"      // class definition
#include "Teuchos_ParameterList.hpp"      // data member
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::Homotopy::Group::Group(
		 Teuchos::ParameterList& locaSublist,
		 const Teuchos::RefCountPtr<LOCA::Homotopy::AbstractGroup>& g,
		 double scalarRandom,
		 double scalarInitialGuess) :
  grpPtr(g),
  gVecPtr(g->getX().clone(NOX::ShapeCopy)),
  randomVecPtr(gVecPtr->clone(NOX::ShapeCopy)),
  newtonVecPtr(),
  gradVecPtr(),
  paramVec(grpPtr->getParams()),
  conParam(0.0),
  conParamID(-1),
  conParamLabel("Homotopy Continuation Parameter"),
  augmentJacForHomotopyNotImplemented(false)
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

LOCA::Homotopy::Group::Group(
		 Teuchos::ParameterList& locaSublist,
		 const Teuchos::RefCountPtr<LOCA::Homotopy::AbstractGroup>& g,
		 const NOX::Abstract::Vector& randomVector) :
  grpPtr(g),
  gVecPtr(g->getX().clone(NOX::ShapeCopy)),
  randomVecPtr(gVecPtr->clone(NOX::ShapeCopy)),
  newtonVecPtr(),
  gradVecPtr(),
  paramVec(grpPtr->getParams()),
  conParam(0.0),
  conParamID(-1),
  conParamLabel("Homotopy Continuation Parameter"),
  augmentJacForHomotopyNotImplemented(false)
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
  grpPtr(Teuchos::rcp_dynamic_cast<LOCA::Homotopy::AbstractGroup>(source.grpPtr->clone(type))),
  gVecPtr(source.gVecPtr->clone(type)),
  randomVecPtr(source.randomVecPtr->clone(NOX::DeepCopy)), //deep copy to keep the set of equations consistent
  newtonVecPtr(),
  gradVecPtr(),
  paramVec(source.paramVec),
  conParam(source.conParam),
  conParamID(source.conParamID),
  conParamLabel(source.conParamLabel),
  augmentJacForHomotopyNotImplemented(source.augmentJacForHomotopyNotImplemented)
{
  if (source.newtonVecPtr != Teuchos::null)
    newtonVecPtr = source.newtonVecPtr->clone(type);

  if (source.gradVecPtr != Teuchos::null)
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
    LOCA::ErrorCheck::throwError("LOCA::Homotopy::Group::Group(copy ctor)",
				 "CopyType is invalid!");
  }

}


LOCA::Homotopy::Group::~Group() 
{
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
    if (newtonVecPtr != Teuchos::null)
      *newtonVecPtr = *(source.newtonVecPtr);
    if (gradVecPtr != Teuchos::null)
      *gradVecPtr = *(source.gradVecPtr);
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;
  }

  return *this;
}

Teuchos::RefCountPtr<NOX::Abstract::Group>
LOCA::Homotopy::Group::clone(NOX::CopyType type) const 
{
  return Teuchos::rcp(new LOCA::Homotopy::Group(*this, type));
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
  string callingFunction = 
    "LOCA::Homotopy::Group::computeDfDp()";
  NOX::Abstract::Group::ReturnType finalStatus;

  // Compute f
  finalStatus = grpPtr->computeF();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // g = conParam * f(x) + ((1.0 - conParam) * (x - randomVec))
  // dg/dp = f(x) - (x - randomVec) where p = conParam
  result = grpPtr->getF();
  result.update(-1.0, grpPtr->getX(), 1.0, *randomVecPtr, 1.0); 

  return finalStatus;
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

  string callingFunction = 
    "LOCA::Homotopy::Group::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus;

  finalStatus = grpPtr->computeF();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  
  // g = conParam * f(x) + ((1.0 - conParam) * (x - randomVec))
  *gVecPtr = grpPtr->getX();
  gVecPtr->update(-1.0, *randomVecPtr, 1.0); 
  gVecPtr->scale( 1.0 - conParam);
  gVecPtr->update(conParam, grpPtr->getF(), 1.0);
 

  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Homotopy::Group::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus;

  finalStatus = grpPtr->computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // Augment the group's Jacobian for homotopy
  NOX::Abstract::Group::ReturnType augHomTest = 
    grpPtr->augmentJacobianForHomotopy(conParam);
  // If it is not implemented, augment the Jacobian during the 
  // applyJacobian() call.
  if (augHomTest == NOX::Abstract::Group::NotDefined)
    augmentJacForHomotopyNotImplemented = true;

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Homotopy::Group::computeGradient()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  finalStatus = computeF();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  
  status = computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  status = applyJacobianTranspose(*gVecPtr, *gradVecPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  return finalStatus;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeNewton(Teuchos::ParameterList& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Homotopy::Group::computeNewton()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  
  if (newtonVecPtr == Teuchos::null)
    newtonVecPtr = gVecPtr->clone(NOX::ShapeCopy);
  
  finalStatus = computeF();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  
  status = computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  status = applyJacobianInverse(params, *gVecPtr, *newtonVecPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  newtonVecPtr->scale(-1.0);
  
  isValidNewton = true;
  
  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobian(const NOX::Abstract::Vector& input,
				   NOX::Abstract::Vector& result) const 
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;
  
  string callingFunction = 
    "LOCA::Homotopy::Group::applyJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus;

  finalStatus = grpPtr->applyJacobian(input, result);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // If the Jacobian is not augmented for homotopy (i.e. using MFNK)
  // then lets augment it here.
  if (augmentJacForHomotopyNotImplemented) {
    double value = 1.0 - conParam;
    result.update(value, input, 1.0);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianTranspose(
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;

  string callingFunction = 
    "LOCA::Homotopy::Group::applyJacobianTranspose()";
  NOX::Abstract::Group::ReturnType finalStatus; 

  finalStatus = grpPtr->applyJacobianTranspose(input, result);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianInverse(Teuchos::ParameterList& params,
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Homotopy::Group::applyJacobianInverse()";
  NOX::Abstract::Group::ReturnType finalStatus; 

  finalStatus = grpPtr->applyJacobianInverse(params, input, result);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  return finalStatus;
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
  if (gradVecPtr == Teuchos::null) {
    LOCA::ErrorCheck::throwError("LOCA::Homotopy::Group::getGradient", 
				 "gradVecPtr is NULL!");
  }
  return *gradVecPtr;
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getNewton() const 
{
  if (newtonVecPtr == Teuchos::null) {
    LOCA::ErrorCheck::throwError("LOCA::Homotopy::Group::getNewton", 
				 "newtonVecPtr is NULL!");
  }
  return *newtonVecPtr;
}

void
LOCA::Homotopy::Group::printSolution(const double conParm) const 
{
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "LOCA::Homotopy::Group::printSolution\n";
    cout << "Nothing to print at this time!";
  }
  return;
}

void
LOCA::Homotopy::Group::printSolution(const NOX::Abstract::Vector& x_,
				     const double conParm) const 
{
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "LOCA::Homotopy::Group::printSolution\n";
    cout << "Nothing to print at this time!";
  }
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
LOCA::Homotopy::Group::setStepperParameters(Teuchos::ParameterList& params)
{
  
  // Create the stepper sublist
  Teuchos::ParameterList& stepperList = params.sublist("Stepper");
  stepperList.set("Continuation Method", "Natural");
  stepperList.set("Continuation Parameter", conParamLabel);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  stepperList.set("Min Value", -1.0);
  stepperList.set("Max Steps", 50);

  // Create predictor sublist
  Teuchos::ParameterList& predictorList = params.sublist("Predictor");
  predictorList.set("Method", "Constant");

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = params.sublist("Step Size");
  //stepSizeList.set("Method", "Constant");
  stepSizeList.set("Method", "Adaptive");
  stepSizeList.set("Initial Step Size", 0.1);
  stepSizeList.set("Min Step Size", 1.0e-2);
  stepSizeList.set("Max Step Size", 1.0);
  stepSizeList.set("Aggressiveness", 0.5);
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
