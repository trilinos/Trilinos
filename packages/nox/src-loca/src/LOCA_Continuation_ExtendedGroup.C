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

#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Parameter_Vector.H"

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
				 LOCA::Continuation::AbstractGroup& grp, 
				 int paramID,
				 NOX::Parameter::List& linSolverParams,
				 NOX::Parameter::List& params)
  : grpPtr(&grp),
    conParamID(paramID),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(false),
    isValidPredictor(false),
    linearSolverParams(linSolverParams),
    theta(params.getParameter("Initial Scale Factor", 1.0)),
    usedConstantPredictor(true)
{
}

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
			      const LOCA::Continuation::AbstractGroup& grp, 
			      int paramID,
			      NOX::Parameter::List& linSolverParams,
			      NOX::Parameter::List& params)
  : grpPtr(dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone())),
    conParamID(paramID),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(true),
    isValidPredictor(false),
    linearSolverParams(linSolverParams),
    theta(params.getParameter("Initial Scale Factor", 1.0)),
    usedConstantPredictor(true)
{
}

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
				 LOCA::Continuation::AbstractGroup& grp, 
				 string paramID,
				 NOX::Parameter::List& linSolverParams,
				 NOX::Parameter::List& params)
  : grpPtr(&grp),
    conParamID(0),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(false),
    isValidPredictor(false),
    linearSolverParams(linSolverParams),
    theta(params.getParameter("Initial Scale Factor", 1.0)),
    usedConstantPredictor(true)
{
  const ParameterVector& p = grpPtr->getParams();
  conParamID = p.getIndex(paramID);
}

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
			     const LOCA::Continuation::AbstractGroup& grp, 
			     string paramID,
			     NOX::Parameter::List& linSolverParams,
			     NOX::Parameter::List& params)
  : grpPtr(dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone())),
    conParamID(0),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(true),
    isValidPredictor(false),
    linearSolverParams(linSolverParams),
    theta(params.getParameter("Initial Scale Factor", 1.0)),
    usedConstantPredictor(true)
{
  const ParameterVector& p = grpPtr->getParams();
  conParamID = p.getIndex(paramID);
}

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
			     const LOCA::Continuation::ExtendedGroup& source, 
			     NOX::CopyType type)
  : grpPtr(dynamic_cast<LOCA::Continuation::AbstractGroup*>(source.grpPtr->clone())),
    conParamID(source.conParamID),
    predictorVec(source.predictorVec),
    ownsGroup(true),
    isValidPredictor(source.isValidPredictor),
    linearSolverParams(source.linearSolverParams),
    theta(source.theta),
    usedConstantPredictor(source.usedConstantPredictor)
{
}


LOCA::Continuation::ExtendedGroup::~ExtendedGroup() 
{
  if (ownsGroup)
    delete grpPtr;
}

NOX::Abstract::Group&
LOCA::Continuation::ExtendedGroup::operator=(
					 const NOX::Abstract::Group& source) 
{
  return operator=(dynamic_cast<const LOCA::Continuation::ExtendedGroup&>(source));
}

LOCA::Extended::AbstractGroup&
LOCA::Continuation::ExtendedGroup::operator=(
				 const LOCA::Extended::AbstractGroup& source) 
{
  return operator=(dynamic_cast<const LOCA::Continuation::ExtendedGroup&>(source));
}

LOCA::Continuation::ExtendedGroup&
LOCA::Continuation::ExtendedGroup::operator=(
			     const LOCA::Continuation::ExtendedGroup& source)
{

  // Protect against A = A
  if (this != &source) {

    // delete old values
    if (ownsGroup)
      delete grpPtr;
    
    // Copy values
    grpPtr = dynamic_cast<LOCA::Continuation::AbstractGroup*>(source.grpPtr->clone());
    conParamID = source.conParamID;
    predictorVec = source.predictorVec;
    ownsGroup = true;
    isValidPredictor = source.isValidPredictor;
    linearSolverParams = source.linearSolverParams;
    theta = source.theta;
    usedConstantPredictor = source.usedConstantPredictor;
  }

  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ExtendedGroup::computeTangent() {
  
  // Get references to x, parameter components of tangent vector
  NOX::Abstract::Vector& tanX = predictorVec.getXVec();
  double& tanP = predictorVec.getParam();

  // Compute Jacobian
  NOX::Abstract::Group::ReturnType res = computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute derivative of residual w.r.t. parameter
  NOX::Abstract::Vector* dfdpVec = tanX.clone(NOX::ShapeCopy);
  res = grpPtr->computeDfDp(conParamID, *dfdpVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Scale dfdp by -1.0
  dfdpVec->scale(-1.0);
  
  // Solve J*tanX = -df/dp
  res = grpPtr->applyJacobianInverse(linearSolverParams, *dfdpVec, tanX);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Set parameter component equal to 1
  tanP = 1.0;

  // If we have a previous solution vector, set orientation equal to 
  // that of the secant vector
  if (isPrevXVec()) {
    // Compute secant vector xVec-prevXVec
    LOCA::Continuation::ExtendedVector secantVec(getPrevX());
    secantVec.update(1.0, getX(), -1.0);
    
    // Give tangent vector same orientation as secant vector
    if (computeScaledDotProduct(secantVec, predictorVec) < 0.0) 
      predictorVec.scale(-1.0);
   
  }

  isValidPredictor = true;

  delete dfdpVec;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::ExtendedGroup::computeSecant() {

  if (!isPrevXVec())
    return NOX::Abstract::Group::Failed;
  
  predictorVec = getPrevX();
  predictorVec.update(1.0, getX(), -1.0);
  predictorVec.scale(1.0/fabs(predictorVec.getParam()));
  
  // else {
//     predictorVec.init(0.0);
//     predictorVec.getParam() = 1.0;
//     usedConstantPredictor = true;
//   }

  isValidPredictor = true;

  return NOX::Abstract::Group::Ok;
}

const LOCA::Continuation::ExtendedVector&
LOCA::Continuation::ExtendedGroup::getPredictorDirection() const {
  return predictorVec;
}

void
LOCA::Continuation::ExtendedGroup::setPredictorDirection(
				const LOCA::Continuation::ExtendedVector& v) 
{
  predictorVec = v;
  isValidPredictor = true;
}

void
LOCA::Continuation::ExtendedGroup::resetIsValid() {
}

bool
LOCA::Continuation::ExtendedGroup::isPredictorDirection() const {
  return isValidPredictor;
}

void
LOCA::Continuation::ExtendedGroup::setContinuationParameter(double val) {
  grpPtr->setParam(conParamID, val);
}

double
LOCA::Continuation::ExtendedGroup::getContinuationParameter() const {
  return grpPtr->getParam(conParamID);
}

const LOCA::Continuation::AbstractGroup&
LOCA::Continuation::ExtendedGroup::getGroup() const {
  return *grpPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Continuation::ExtendedGroup::getGroup() {
  return *grpPtr;
}

void
LOCA::Continuation::ExtendedGroup::setScaleFactor(double t) {
  theta = t;
}

double
LOCA::Continuation::ExtendedGroup::getScaleFactor() const {
  return theta;
}

double
LOCA::Continuation::ExtendedGroup::getStepSizeScaleFactor() const {
  return 1.0;
}

double
LOCA::Continuation::ExtendedGroup::computeScaledDotProduct(
			   const LOCA::Continuation::ExtendedVector& x,
			   const LOCA::Continuation::ExtendedVector& y) const
{
  double val = grpPtr->computeScaledDotProduct(x.getXVec(), y.getXVec());
  return val + theta*theta*x.getParam()*y.getParam();
}

double
LOCA::Continuation::ExtendedGroup::computeScaledDotProduct(
				     const NOX::Abstract::Vector& x,
				     const NOX::Abstract::Vector& y) const
{
  return computeScaledDotProduct(
		   dynamic_cast<const LOCA::Continuation::ExtendedVector&>(x),
		   dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y));
}

void
LOCA::Continuation::ExtendedGroup::printSolution() const
{
  grpPtr->printSolution(getContinuationParameter());
}

const LOCA::Continuation::AbstractGroup&
LOCA::Continuation::ExtendedGroup::getUnderlyingGroup() const
{
  return *grpPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Continuation::ExtendedGroup::getUnderlyingGroup()
{
  return *grpPtr;
}
