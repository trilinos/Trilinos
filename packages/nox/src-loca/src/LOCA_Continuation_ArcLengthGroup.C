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

#include "LOCA_Continuation_ArcLengthGroup.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"

LOCA::Continuation::ArcLengthGroup::ArcLengthGroup(
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
    arclengthStep(0.0),
    isValidPrevXVec(false),
    stepSizeScaleFactor(1.0),
    isFirstRescale(true)
{
  resetIsValid();
  doArcLengthScaling = params.getParameter("Enable Arc Length Scaling",true); 
  gGoal = params.getParameter("Goal Arc Length Parameter Contribution", 0.5);
  gMax = params.getParameter("Max Arc Length Parameter Contribution", 0.0);
  thetaMin = params.getParameter("Min Scale Factor", 1.0e-3);
}

LOCA::Continuation::ArcLengthGroup::ArcLengthGroup(
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
    arclengthStep(0.0),
    isValidPrevXVec(false),
    stepSizeScaleFactor(1.0),
    isFirstRescale(true)
{
  resetIsValid();
  doArcLengthScaling = params.getParameter("Enable Arc Length Scaling",true);
  gGoal = params.getParameter("Goal Arc Length Parameter Contribution", 0.5);
  gMax = params.getParameter("Max Arc Length Parameter Contribution", 0.0);
  thetaMin = params.getParameter("Min Scale Factor", 1.0e-3);
}

LOCA::Continuation::ArcLengthGroup::ArcLengthGroup(
			     const LOCA::Continuation::ArcLengthGroup& source,
			     NOX::CopyType type)
  :  LOCA::Continuation::ExtendedGroup(source, type), 
    xVec(source.xVec, type),
    fVec(source.fVec, type),
    newtonVec(source.newtonVec, type),
    gradientVec(source.gradientVec, type),
    prevXVec(source.prevXVec, type),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    arclengthStep(source.arclengthStep),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    isValidPrevXVec(source.isValidPrevXVec),
    doArcLengthScaling(source.doArcLengthScaling),
    gGoal(source.gGoal),
    gMax(source.gMax),
    thetaMin(source.thetaMin),
    stepSizeScaleFactor(source.stepSizeScaleFactor),
    isFirstRescale(source.isFirstRescale)
{
}


LOCA::Continuation::ArcLengthGroup::~ArcLengthGroup() 
{
  delete derivResidualParamPtr;
}

LOCA::Continuation::ExtendedGroup&
LOCA::Continuation::ArcLengthGroup::operator=(
			      const LOCA::Continuation::ExtendedGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Continuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Continuation::ArcLengthGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::Continuation::ArcLengthGroup&>(source);
}

LOCA::Continuation::ArcLengthGroup&
LOCA::Continuation::ArcLengthGroup::operator=(
			    const LOCA::Continuation::ArcLengthGroup& source) 
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
    arclengthStep = source.arclengthStep;
    *derivResidualParamPtr = *source.derivResidualParamPtr;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidPrevXVec = source.isValidPrevXVec;
    doArcLengthScaling = source.doArcLengthScaling;
    gGoal = source.gGoal;
    gMax = source.gMax;
    thetaMin = source.thetaMin;
    stepSizeScaleFactor = source.stepSizeScaleFactor;
    isFirstRescale = source.isFirstRescale;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Continuation::ArcLengthGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Continuation::ArcLengthGroup(*this);
}

void
LOCA::Continuation::ArcLengthGroup::setStepSize(double deltaS) 
{
  arclengthStep = deltaS;
  
  // Arclength Step appears on RHS but not in Jacobian
  isValidF = false;
  isValidNewton = false;
  isValidGradient = false;
}

double
LOCA::Continuation::ArcLengthGroup::getStepSize() const 
{
  return arclengthStep;
}

void 
LOCA::Continuation::ArcLengthGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y) );
}

void
LOCA::Continuation::ArcLengthGroup::setX(
				 const LOCA::Continuation::ExtendedVector& y) 
{
  grpPtr->setX( y.getXVec() );
  grpPtr->setParam(conParamID, y.getParam());
  xVec = y;

  resetIsValid();
}

void
LOCA::Continuation::ArcLengthGroup::setPrevX(const NOX::Abstract::Vector& y) 
{
  setPrevX( dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y) );
}

void
LOCA::Continuation::ArcLengthGroup::setPrevX(
				 const LOCA::Continuation::ExtendedVector& y) 
{
  prevXVec = y;

  isValidF = false; // Previous vector part of residual equation
  isValidNewton = false;
  isValidGradient = false;
  isValidPrevXVec = true;
}

const LOCA::Continuation::ExtendedVector&
LOCA::Continuation::ArcLengthGroup::getPrevX() const 
{
  return prevXVec;
}

bool
LOCA::Continuation::ArcLengthGroup::isPrevXVec() const
{
  return isValidPrevXVec;
}

void
LOCA::Continuation::ArcLengthGroup::computeX(const NOX::Abstract::Group& g, 
					     const NOX::Abstract::Vector& d,
					     double step) 
{
  computeX( dynamic_cast<const LOCA::Continuation::ArcLengthGroup&>(g),
	    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(d),
	    step);
}

void
LOCA::Continuation::ArcLengthGroup::computeX(
			       const LOCA::Continuation::ArcLengthGroup& g, 
			       const LOCA::Continuation::ExtendedVector& d,
			       double step) 
{
  grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  xVec.update(1.0, g.getX(), step, d, 0.0);
  grpPtr->setParam(conParamID, xVec.getParam());

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  string callingFunction = "LOCA::Continuation::ArcLengthGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // make sure predictor is valid
  if (!isPredictorDirection()) {
    LOCA::ErrorCheck::throwError(callingFunction, 
				 "Called with invalid predictor vector.");
  }

  // Compute underlying F
  if (!grpPtr->isF()) {
    finalStatus = grpPtr->computeF();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  fVec.getXVec() = grpPtr->getF();
  
  // Compute xVec-prevXVec;
  LOCA::Continuation::ExtendedVector *tmpVec =
    dynamic_cast<LOCA::Continuation::ExtendedVector*>(xVec.clone(NOX::DeepCopy));
  tmpVec->update(-1.0, prevXVec, 1.0);
  
  fVec.getParam() =  
    computeScaledDotProduct(predictorVec, *tmpVec) - arclengthStep * computeScaledDotProduct(predictorVec, predictorVec);

  delete tmpVec;
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Continuation::ArcLengthGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  status = grpPtr->computeDfDp(conParamID, *derivResidualParamPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

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
LOCA::Continuation::ArcLengthGroup::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Continuation::ArcLengthGroup::computeGradient()";
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
  gradientVec.getParam() = derivResidualParamPtr->dot(fVec.getXVec());
  gradientVec.update(fVec.getParam(), predictorVec, 1.0);

  isValidGradient = true;

  return finalStatus;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeNewton(
					       NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Continuation::ArcLengthGroup::computeNewton()";
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

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    double r = computeScaledDotProduct(newtonVec, predictorVec);
    cout << "\n\tScaled component of Newton vector in direction of "
	 << "predictor:  " << r << endl;
  }

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::applyJacobian(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Continuation::ArcLengthGroup::applyJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // make sure predictor is valid
  if (!isPredictorDirection()) {
    LOCA::ErrorCheck::throwError(callingFunction, 
				 "Called with invalid predictor vector.");
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

  // compute J*x
  status = grpPtr->applyJacobian(input_x, result_x);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute J*x + p*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, 1.0);

  // compute dx/ds x + dp/ds p
  result_param = computeScaledDotProduct(predictorVec, c_input);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input, 
					  NOX::Abstract::Vector& result) const 
{

  string callingFunction = 
    "LOCA::Continuation::ArcLengthGroup::applyJacobianTranspose()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // make sure predictor is valid
  if (!isPredictorDirection()) {
    LOCA::ErrorCheck::throwError(callingFunction, 
				 "Called with invalid predictor vector.");
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
  result_param = derivResidualParamPtr->dot(input_x);

  // compute [J^T x; df/dp^T x] + input_p * [dx/ds; dp/ds]
  c_result.update(input_param, predictorVec, 1.0);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::applyJacobianInverse(
					  NOX::Parameter::List& params, 
					  const NOX::Abstract::Vector& input, 
					  NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Continuation::ArcLengthGroup::applyJacobianInverse()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // make sure predictor is valid
  if (!isPredictorDirection()) {
    LOCA::ErrorCheck::throwError(callingFunction, 
				 "Called with invalid predictor vector.");
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

  const NOX::Abstract::Vector** rhs = new const NOX::Abstract::Vector*[2];
  NOX::Abstract::Vector** lhs = new NOX::Abstract::Vector*[2];
  rhs[0] = &input_x;
  rhs[1] = derivResidualParamPtr;
  lhs[0] = input_x.clone(NOX::ShapeCopy);
  lhs[1] = input_x.clone(NOX::ShapeCopy);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Solve J*lhs = rhs
  status = grpPtr->applyJacobianInverseMulti(params, rhs, lhs, 2);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Get x, param components of predictor vector
  const NOX::Abstract::Vector& tanX = predictorVec.getXVec();
  double tanP = predictorVec.getParam();

  // Compute result_param
  result_param =
    (grpPtr->computeScaledDotProduct(tanX, *lhs[0]) - input_param) / 
    (grpPtr->computeScaledDotProduct(tanX, *lhs[1]) - theta*theta*tanP);

  // Compute result_x = a + result_param*b 
  result_x.update(1.0, *lhs[0], -result_param, *lhs[1], 0.0);

  // Clean up memory
  delete lhs[0];
  delete lhs[1];
  delete [] rhs;
  delete [] lhs;

  return finalStatus;
}

bool
LOCA::Continuation::ArcLengthGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Continuation::ArcLengthGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Continuation::ArcLengthGroup::isGradient() const 
{
  return isValidGradient;
}

bool
LOCA::Continuation::ArcLengthGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Continuation::ArcLengthGroup::getX() const 
{
  return xVec;
}

const NOX::Abstract::Vector&
LOCA::Continuation::ArcLengthGroup::getF() const 
{
  return fVec;
}

double
LOCA::Continuation::ArcLengthGroup::getNormF() const 
{
  return fVec.norm();
}

const NOX::Abstract::Vector&
LOCA::Continuation::ArcLengthGroup::getGradient() const 
{
  return gradientVec;
}

const NOX::Abstract::Vector&
LOCA::Continuation::ArcLengthGroup::getNewton() const 
{
  return newtonVec;
}

double
LOCA::Continuation::ArcLengthGroup::getNormNewtonSolveResidual() const 
{
  string callingFunction = 
    "LOCA::Continuation::ArcLengthGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::Continuation::ExtendedVector residual = fVec;
  
  finalStatus = applyJacobian(newtonVec, residual);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  residual = residual.update(1.0, fVec, 1.0);
  return residual.norm();
}

void
LOCA::Continuation::ArcLengthGroup::setContinuationParameter(double val) 
{
  LOCA::Continuation::ExtendedGroup::setContinuationParameter(val);
  xVec.getParam() = val;
  resetIsValid();
}

void
LOCA::Continuation::ArcLengthGroup::resetIsValid() {
  LOCA::Continuation::ExtendedGroup::resetIsValid();
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
}

void
LOCA::Continuation::ArcLengthGroup::scalePredictor(
				       LOCA::Continuation::ExtendedVector& v) {

  if (doArcLengthScaling) {

    // Estimate dpds
    double dpdsOld = 
      1.0/sqrt(computeScaledDotProduct(v, v));

    if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
      cout << endl 
	   << "\t" << LOCA::Utils::fill(64, '+') << endl 
	   << "\t" << "Arc-length scaling calculation:" << endl
	   << "\t" << "Parameter component of predictor before rescaling = " 
	   << LOCA::Utils::sci(dpdsOld) << endl
	   << "\t" << "Scale factor from previous step                   = "
	   << LOCA::Utils::sci(theta) << endl
	   << "\t" << "Parameter contribution to arc-length equation     = "
	   << LOCA::Utils::sci(theta*dpdsOld) << endl;
    }

    // Recompute scale factor
    recalculateScaleFactor(dpdsOld);

    // Calculate new dpds using new scale factor
    double dpdsNew = 
      1.0/sqrt(computeScaledDotProduct(v, v));

    if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
      cout << endl 
	   << "\t" << "Parameter component of predictor after rescaling  = " 
	   << LOCA::Utils::sci(dpdsNew) << endl
	   << "\t" << "New scale factor (theta)                          = "
	   << LOCA::Utils::sci(theta) << endl
	   << "\t" << "Parameter contribution to arc-length equation     = "
	   << LOCA::Utils::sci(theta*dpdsNew) << endl
	   << "\t" << LOCA::Utils::fill(64, '+') << endl;
    }

    // Rescale predictor vector
    v.scale(dpdsNew);

    // Adjust step size scaling factor to reflect changes in 
    // arc-length parameterization
    // The first time we rescale (first continuation step) we use a different
    // step size scale factor so that dpds*deltaS = step size provided by user
    if (isFirstRescale) {
      stepSizeScaleFactor = 1.0/dpdsNew;
      isFirstRescale = false;
    }
    else
      stepSizeScaleFactor = dpdsOld/dpdsNew;
  }

  return;
}

void
LOCA::Continuation::ArcLengthGroup::recalculateScaleFactor(double dpds) {
  
  double thetaOld = getScaleFactor();
  double g = dpds*thetaOld;

  if (g > gMax) {
    double thetaNew;
    
    thetaNew = gGoal/dpds * sqrt( (1.0 - g*g) / (1.0 - gGoal*gGoal) ); 

    if (thetaNew < thetaMin)
      thetaNew = thetaMin;

    setScaleFactor(thetaNew);
  }

}

double
LOCA::Continuation::ArcLengthGroup::getStepSizeScaleFactor() const {
  return stepSizeScaleFactor;
}


