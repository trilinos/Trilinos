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

#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Parameter_Vector.H"

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
				 LOCA::Continuation::AbstractGroup& grp, 
				 int paramID,
				 NOX::Parameter::List& params)
  : grpPtr(&grp),
    conParamID(paramID),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(false),
    isValidPredictor(false),
    theta(params.getParameter("Initial Scale Factor", 1.0))
{
}

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
			      const LOCA::Continuation::AbstractGroup& grp, 
			      int paramID,
			      NOX::Parameter::List& params)
  : grpPtr(dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone())),
    conParamID(paramID),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(true),
    isValidPredictor(false),
    theta(params.getParameter("Initial Scale Factor", 1.0))
{
}

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
				 LOCA::Continuation::AbstractGroup& grp, 
				 string paramID,
				 NOX::Parameter::List& params)
  : grpPtr(&grp),
    conParamID(0),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(false),
    isValidPredictor(false),
    theta(params.getParameter("Initial Scale Factor", 1.0))
{
  const ParameterVector& p = grpPtr->getParams();
  conParamID = p.getIndex(paramID);
}

LOCA::Continuation::ExtendedGroup::ExtendedGroup(
			     const LOCA::Continuation::AbstractGroup& grp, 
			     string paramID,
			     NOX::Parameter::List& params)
  : grpPtr(dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone())),
    conParamID(0),
    predictorVec(grp.getX(), 0.0),
    ownsGroup(true),
    isValidPredictor(false),
    theta(params.getParameter("Initial Scale Factor", 1.0))
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
    theta(source.theta)
{
}


LOCA::Continuation::ExtendedGroup::~ExtendedGroup() 
{
  if (ownsGroup) {
    delete grpPtr;
  }
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
    *grpPtr = *(source.grpPtr);
    conParamID = source.conParamID;
    predictorVec = source.predictorVec;
    isValidPredictor = source.isValidPredictor;
    theta = source.theta;
  }

  return *this;
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

int
LOCA::Continuation::ExtendedGroup::getContinuationParameterID() const {
  return conParamID;
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
				     const NOX::Abstract::Vector& x,
				     const NOX::Abstract::Vector& y) const
{
  return computeScaledDotProduct(
		   dynamic_cast<const LOCA::Continuation::ExtendedVector&>(x),
		   dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y));
}

double
LOCA::Continuation::ExtendedGroup::computeScaledDotProduct(
			   const LOCA::Continuation::ExtendedVector& x,
			   const LOCA::Continuation::ExtendedVector& y) const
{
  double val = grpPtr->computeScaledDotProduct(x.getXVec(), y.getXVec());
  return val + theta*theta*x.getParam()*y.getParam();
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
