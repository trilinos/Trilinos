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

#include "LOCA_Continuation_ArcLengthGroup.H"
#include "LOCA_Utils.H"

LOCA::Continuation::ArcLengthGroup::ArcLengthGroup(
				 LOCA::Continuation::AbstractGroup& g,
				 int paramID,
				 NOX::Parameter::List& linSolverParams,
				 NOX::Parameter::List& params)
  : LOCA::Continuation::ExtendedGroup(g, paramID, linSolverParams, params), 
    xVec(g.getX(), g.getParam(paramID)),
    fVec(g.getX(), 0.0),
    newtonVec(g.getX(), 0.0),
    prevXVec(g.getX(), g.getParam(paramID)),
    derivResidualParamPtr(g.getX().clone(NOX::ShapeCopy)),
    arclengthStep(0.0),
    isValidPrevXVec(false),
    stepSizeScaleFactor(1.0)
{
  resetIsValid();
  gGoal = params.getParameter("Goal g", 0.5);
  gMax = params.getParameter("Max g", 0.0);
  thetaMin = params.getParameter("Min Scale Factor", 1.0e-3);
}

LOCA::Continuation::ArcLengthGroup::ArcLengthGroup(
				  LOCA::Continuation::AbstractGroup& g,
				  string paramID,
				  NOX::Parameter::List& linSolverParams,
				  NOX::Parameter::List& params)
  : LOCA::Continuation::ExtendedGroup(g, paramID, linSolverParams, params), 
    xVec(g.getX(), g.getParam(paramID)),
    fVec(g.getX(), 0.0),
    newtonVec(g.getX(), 0.0),
    prevXVec(g.getX(), g.getParam(paramID)),
    derivResidualParamPtr(g.getX().clone(NOX::ShapeCopy)),
    arclengthStep(0.0),
    isValidPrevXVec(false),
    stepSizeScaleFactor(1.0)
{
  resetIsValid();
  gGoal = params.getParameter("Goal g", 0.5);
  gMax = params.getParameter("Max g", 0.0);
  thetaMin = params.getParameter("Min Scale Factor", 1.0e-3);
}

LOCA::Continuation::ArcLengthGroup::ArcLengthGroup(const LOCA::Continuation::ArcLengthGroup& source, NOX::CopyType type)
  :  LOCA::Continuation::ExtendedGroup(source, type), 
    xVec(source.xVec, type),
    fVec(source.fVec, type),
    newtonVec(source.newtonVec, type),
    prevXVec(source.prevXVec, type),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    arclengthStep(source.arclengthStep),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    isValidPrevXVec(source.isValidPrevXVec),
    gGoal(source.gGoal),
    gMax(source.gMax),
    thetaMin(source.thetaMin),
    stepSizeScaleFactor(source.stepSizeScaleFactor)
{
}


LOCA::Continuation::ArcLengthGroup::~ArcLengthGroup() 
{
  delete derivResidualParamPtr;
}

LOCA::Continuation::ExtendedGroup&
LOCA::Continuation::ArcLengthGroup::operator=(const LOCA::Continuation::ExtendedGroup& source)
{
  return *this = dynamic_cast<const LOCA::Continuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Continuation::ArcLengthGroup::operator=(const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const LOCA::Continuation::ArcLengthGroup&>(source);
}

LOCA::Continuation::ArcLengthGroup&
LOCA::Continuation::ArcLengthGroup::operator=(const LOCA::Continuation::ArcLengthGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    LOCA::Continuation::ExtendedGroup::operator=(source);

    // Delete old values
    delete derivResidualParamPtr;

    // Copy values
    xVec = source.xVec;
    fVec = source.fVec;
    newtonVec = source.newtonVec;
    prevXVec = source.prevXVec;
    arclengthStep = source.arclengthStep;
    derivResidualParamPtr = source.derivResidualParamPtr->clone(NOX::DeepCopy);
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidPrevXVec = source.isValidPrevXVec;
    gGoal = source.gGoal;
    gMax = source.gMax;
    thetaMin = source.thetaMin;
    stepSizeScaleFactor = source.stepSizeScaleFactor;
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
LOCA::Continuation::ArcLengthGroup::setX(const LOCA::Continuation::ExtendedVector& y) 
{
  LOCA::Continuation::ExtendedGroup::grpPtr->setX( y.getXVec() );
  LOCA::Continuation::ExtendedGroup::grpPtr->setParam(LOCA::Continuation::ExtendedGroup::conParamID, y.getParam());
  xVec = y;

  resetIsValid();
}

void
LOCA::Continuation::ArcLengthGroup::setPrevX(const NOX::Abstract::Vector& y) 
{
  setPrevX( dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y) );
}

void
LOCA::Continuation::ArcLengthGroup::setPrevX(const LOCA::Continuation::ExtendedVector& y) 
{
  prevXVec = y;

  isValidF = false; // Previous vector part of residual equation
  isValidNewton = false;
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
LOCA::Continuation::ArcLengthGroup::computeX(const LOCA::Continuation::ArcLengthGroup& g, 
					     const LOCA::Continuation::ExtendedVector& d,
					     double step) 
{
  LOCA::Continuation::ExtendedGroup::grpPtr->computeX(*(g.LOCA::Continuation::ExtendedGroup::grpPtr), d.getXVec(), step);
  xVec.update(1.0, g.getX(), step, d, 0.0);
  LOCA::Continuation::ExtendedGroup::grpPtr->setParam(LOCA::Continuation::ExtendedGroup::conParamID, xVec.getParam());

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  res = LOCA::Continuation::ExtendedGroup::grpPtr->computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  fVec.getXVec() = LOCA::Continuation::ExtendedGroup::grpPtr->getF();
  
  // Compute xVec-prevXVec;
  LOCA::Continuation::ExtendedVector *tmpVec =
    dynamic_cast<LOCA::Continuation::ExtendedVector*>(xVec.clone(NOX::DeepCopy));
  tmpVec->update(-1.0, prevXVec, 1.0);

  if (!isPredictorDirection()) {
    cerr << "LOCA::Continuation::ArcLengthGroup::computeF() called " 
	 << "with invalid predictor vector." << endl;
    return NOX::Abstract::Group::Failed;
  }

  fVec.getParam() =  
    computeScaledDotProduct(predictorVec, *tmpVec) - arclengthStep;

  delete tmpVec;
  
  isValidF = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  res = LOCA::Continuation::ExtendedGroup::grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = LOCA::Continuation::ExtendedGroup::grpPtr->computeDfDp(LOCA::Continuation::ExtendedGroup::conParamID, *derivResidualParamPtr);

  if (res != NOX::Abstract::Group::Ok)
    return res;

  isValidJacobian = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeNewton(NOX::Parameter::List& params) 
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
LOCA::Continuation::ArcLengthGroup::computeTangent() {
  
  // Compute tangent vector using base class which assumes |dp/ds| = 1
  NOX::Abstract::Group::ReturnType res =
    LOCA::Continuation::ExtendedGroup::computeTangent();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  scalePredictor();

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::computeSecant() {
  
  // Compute secant vector using base class
  NOX::Abstract::Group::ReturnType res =
    LOCA::Continuation::ExtendedGroup::computeSecant();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  if (isPrevXVec())
    scalePredictor();

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
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

  NOX::Abstract::Group::ReturnType res;

  // compute J*x
  res = LOCA::Continuation::ExtendedGroup::grpPtr->applyJacobian(input_x, result_x);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*x + p*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, 1.0);

  // if tangent vector hasn't been computed, we are stuck since this is
  // a const method
  if (!isPredictorDirection()) {
    cerr << "LOCA::Continuation::ArcLengthGroup::applyJacobian() called " 
	 << "with invalid predictor vector." << endl;
    return NOX::Abstract::Group::Failed;
  }

  // compute dx/ds x + dp/ds p
  result_param = computeScaledDotProduct(predictorVec, c_input);

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ArcLengthGroup::applyJacobianInverse(NOX::Parameter::List& params, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
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
  
  NOX::Abstract::Group::ReturnType res;

  const NOX::Abstract::Vector** rhs = new const NOX::Abstract::Vector*[2];
  NOX::Abstract::Vector** lhs = new NOX::Abstract::Vector*[2];
  rhs[0] = &input_x;
  rhs[1] = derivResidualParamPtr;
  lhs[0] = input_x.clone(NOX::ShapeCopy);
  lhs[1] = input_x.clone(NOX::ShapeCopy);

  // Solve J*lhs = rhs
  res = grpPtr->applyJacobianInverseMulti(params, rhs, lhs, 2);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Get x, param components of predictor vector
  const NOX::Abstract::Vector& tanX =  
    LOCA::Continuation::ExtendedGroup::predictorVec.getXVec();
  double tanP = 
    LOCA::Continuation::ExtendedGroup::predictorVec.getParam();

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

  return res;
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
  return false;
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
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Continuation::ArcLengthGroup::getNewton() const 
{
  return newtonVec;
}

double
LOCA::Continuation::ArcLengthGroup::getNormNewtonSolveResidual() const 
{
  LOCA::Continuation::ExtendedVector residual = fVec;
  
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
}

void
LOCA::Continuation::ArcLengthGroup::scalePredictor() {

  // Estimate dpds
  double dpdsOld = 1.0/sqrt(computeScaledDotProduct(predictorVec, predictorVec));

  if (Utils::doPrint(Utils::StepperDetails)) {
    cout << "LOCA::Continuation::ArcLengthGroup::scalePredictor():  "
         << "dpdsOld = " << dpdsOld << endl;
    cout << "LOCA::Continuation::ArcLengthGroup::scalePredictor():  "
         << "thetaOld = " << theta << endl;
    cout << "LOCA::Continuation::ArcLengthGroup::scalePredictor():  "
         << "gOld = " << theta*dpdsOld << endl;
  }

  // Recompute scale factor
  recalculateScaleFactor(dpdsOld);

  // Calculate new dpds using new scale factor
  double dpdsNew = 1.0/sqrt(computeScaledDotProduct(predictorVec, predictorVec));

  if (Utils::doPrint(Utils::StepperDetails)) {
    cout << "LOCA::Continuation::ArcLengthGroup::scalePredictor():  "
         << "dpdsNew = " << dpdsNew << endl;
    cout << "LOCA::Continuation::ArcLengthGroup::scalePredictor():  "
         << "thetaNew = " << theta << endl;
    cout << "LOCA::Continuation::ArcLengthGroup::scalePredictor():  "
         << "gNew = " << theta*dpdsNew << endl;
  }

  // Rescale predictor vector
  predictorVec.scale(dpdsNew);

  // Adjust step size scaling factor to reflect changes in 
  // arc-length parameterization
  if (LOCA::Continuation::ExtendedGroup::usedConstantPredictor && isPrevXVec()) {
    stepSizeScaleFactor = 1.0/dpdsNew;
    LOCA::Continuation::ExtendedGroup::usedConstantPredictor = false;
  }
  else
    stepSizeScaleFactor = dpdsOld/dpdsNew;

  return;
}

void
LOCA::Continuation::ArcLengthGroup::recalculateScaleFactor(double dpds) {
  
  double thetaOld = getScaleFactor();
  double g = dpds*thetaOld;

  static bool isFirstIt = false;

  // if (g >= 0.9)
//     return;

  if (g > gMax || isFirstIt) {
    double thetaNew;
    
    thetaNew = gGoal/dpds * sqrt( (1.0 - g*g) / (1.0 - gGoal*gGoal) ); 

    if (thetaNew < thetaMin)
      thetaNew = thetaMin;

    setScaleFactor(thetaNew);

    if (isFirstIt)
      isFirstIt = false;
  }

}

double
LOCA::Continuation::ArcLengthGroup::getStepSizeScaleFactor() const {
  return stepSizeScaleFactor;
}


