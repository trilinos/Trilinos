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

#include "LOCA_Continuation_Group.H"
#include "LOCA_Parameter_Vector.H"

LOCA::Continuation::Group::Group(const LOCA::Abstract::Group& grp, 
				 int paramID,
				 const NOX::Parameter::List& linSolverParams,
				 NOX::Parameter::List& params)
  : grpPtr(dynamic_cast<LOCA::Abstract::Group*>(grp.clone(NOX::DeepCopy))),
    conParamID(paramID),
    tangentVec(grp.getX(), 0.0),
    isValidTangent(false),
    linearSolverParams(linSolverParams),
    scaleVec(grp.getScaleVec(), 
	     params.getParameter("Initial Scale Factor", 1.0))
{
}

LOCA::Continuation::Group::Group(const LOCA::Abstract::Group& grp, 
				 string paramID,
				 const NOX::Parameter::List& linSolverParams,
				 NOX::Parameter::List& params)
  : grpPtr(dynamic_cast<LOCA::Abstract::Group*>(grp.clone(NOX::DeepCopy))),
    conParamID(0),
    tangentVec(grp.getX(), 0.0),
    isValidTangent(false),
    linearSolverParams(linSolverParams),
    scaleVec(grp.getScaleVec(), 
	     params.getParameter("Initial Scale Factor", 1.0))
{
  const ParameterVector& p = grpPtr->getParams();
  conParamID = p.getIndex(paramID);
}

LOCA::Continuation::Group::Group(const Group& source, NOX::CopyType type)
  : grpPtr(dynamic_cast<LOCA::Abstract::Group*>(source.grpPtr->clone(type))),
    conParamID(source.conParamID),
    tangentVec(source.tangentVec),
    isValidTangent(source.isValidTangent),
    linearSolverParams(source.linearSolverParams),
    scaleVec(source.scaleVec)
{
}


LOCA::Continuation::Group::~Group() 
{
  delete grpPtr;
}

LOCA::Continuation::Group&
LOCA::Continuation::Group::operator=(const LOCA::Continuation::Group& source)
{

  // Protect against A = A
  if (this != &source) {
    NOX::CopyType type = NOX::DeepCopy;

    // Delete old values
    delete grpPtr;

    // Copy values
    grpPtr = dynamic_cast<LOCA::Abstract::Group*>(source.grpPtr->clone(type));
    conParamID = source.conParamID;
    tangentVec = source.tangentVec;
    isValidTangent = source.isValidTangent;
    linearSolverParams = source.linearSolverParams;
    scaleVec = source.scaleVec;
  }

  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::Group::computeTangent(NOX::Parameter::List& params) {
  
  // Get references to x, parameter components of tangent vector
  NOX::Abstract::Vector& tanX = tangentVec.getXVec();
  double& tanP = tangentVec.getParam();

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
  res = grpPtr->applyJacobianInverse(params, *dfdpVec, tanX);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Set parameter component equal to 1
  tanP = 1.0;

  // If we have a previous solution vector, set orientation equal to 
  // that of the secant vector
  if (isPrevXVec()) {
    // Compute secant vector xVec-prevXVec
    LOCA::Continuation::Vector secantVec(getPrevX());
    secantVec.update(1.0, getX(), -1.0);
    
    // Give tangent vector same orientation as secant vector
    if (scaledDotProduct(secantVec, tangentVec) < 0.0) 
      tangentVec.scale(-1.0);
   
  }

  isValidTangent = true;

  delete dfdpVec;

  return res;
}

const LOCA::Continuation::Vector&
LOCA::Continuation::Group::getTangent() const {
  return tangentVec;
}

void
LOCA::Continuation::Group::resetIsValid() {
}

bool
LOCA::Continuation::Group::isTangent() const {
  return isValidTangent;
}

void
LOCA::Continuation::Group::setContinuationParameter(double val) {
  grpPtr->setParam(conParamID, val);
}

double
LOCA::Continuation::Group::getContinuationParameter() const {
  return grpPtr->getParam(conParamID);
}

const LOCA::Abstract::Group&
LOCA::Continuation::Group::getGroup() const {
  return *grpPtr;
}

void
LOCA::Continuation::Group::setScaleFactor(double theta) {
  scaleVec.getParam() = theta;
}

double
LOCA::Continuation::Group::getScaleFactor() const {
  return scaleVec.getParam();
}

double
LOCA::Continuation::Group::scaledDotProduct(const LOCA::Continuation::Vector& x,
					    const LOCA::Continuation::Vector& y) const
{

  return LOCA::Continuation::scaledDotProduct(x, y, scaleVec);
}

double
LOCA::Continuation::scaledDotProduct(const NOX::Abstract::Vector& x,
				     const NOX::Abstract::Vector& y,
				     const NOX::Abstract::Vector& s) 
{

  double d;

  // Copy x
  NOX::Abstract::Vector* sxPtr = x.clone();

  // Scale x twice
  sxPtr->scale(s);
  sxPtr->scale(s);

  // Compute dot product
  d = sxPtr->dot(y);

  // Delete copy
  delete sxPtr;

  return d;
}
