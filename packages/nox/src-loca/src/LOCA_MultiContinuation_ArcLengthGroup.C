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

#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "NOX_Parameter_List.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
				 LOCA::MultiContinuation::AbstractGroup& g, 
				 const vector<int>& paramIDs,
				 NOX::Parameter::List& params)
  : LOCA::MultiContinuation::ExtendedGroup(g, paramIDs, params),
    constraints(numParams, 1),
    isValidConstraints(false),
    isValidConstraintDerivatives(false),
    theta(numParams, 1.0),
    isFirstRescale(true)
{
  double theta0 = params.getParameter("Initial Scale Factor", 1.0);
  doArcLengthScaling = params.getParameter("Enable Arc Length Scaling",true); 
  gGoal = params.getParameter("Goal Arc Length Parameter Contribution", 0.5);
  gMax = params.getParameter("Max Arc Length Parameter Contribution", 0.0);
  thetaMin = params.getParameter("Min Scale Factor", 1.0e-3);
  
  for (int i=0; i<numParams; i++)
    theta[i] = theta0;
}

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
				 LOCA::MultiContinuation::AbstractGroup& g, 
				 const string& paramID,
				 NOX::Parameter::List& params)
  : LOCA::MultiContinuation::ExtendedGroup(g, paramID, params),
    constraints(numParams, 1),
    isValidConstraints(false),
    isValidConstraintDerivatives(false),
    theta(numParams, 1.0),
    isFirstRescale(true)
{
  double theta0 = params.getParameter("Initial Scale Factor", 1.0);
  doArcLengthScaling = params.getParameter("Enable Arc Length Scaling",true); 
  gGoal = params.getParameter("Goal Arc Length Parameter Contribution", 0.5);
  gMax = params.getParameter("Max Arc Length Parameter Contribution", 0.0);
  thetaMin = params.getParameter("Min Scale Factor", 1.0e-3);
  
  for (int i=0; i<numParams; i++)
    theta[i] = theta0;
}

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
			 const LOCA::MultiContinuation::ArcLengthGroup& source,
			 NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type),
    constraints(source.constraints),
    isValidConstraints(source.isValidConstraints),
    isValidConstraintDerivatives(source.isValidConstraintDerivatives),
    theta(source.theta),
    doArcLengthScaling(source.doArcLengthScaling),
    gGoal(source.gGoal),
    gMax(source.gMax),
    thetaMin(source.thetaMin),
    isFirstRescale(source.isFirstRescale)
{
}


LOCA::MultiContinuation::ArcLengthGroup::~ArcLengthGroup() 
{
}

LOCA::MultiContinuation::ArcLengthGroup&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
		       const LOCA::MultiContinuation::ArcLengthGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    LOCA::MultiContinuation::ExtendedGroup::operator=(source);
    constraints = source.constraints;
    isValidConstraints = source.isValidConstraints;
    isValidConstraintDerivatives = source.isValidConstraintDerivatives;
    theta = source.theta;
  }

  return *this;
}

LOCA::MultiContinuation::ExtendedGroup&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
		        const LOCA::MultiContinuation::ExtendedGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthGroup::computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::NultiContinuation::ArcLengthGroup::computeConstraintss()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // compute predictor
  if (!isValidConstraintDerivatives) {
    finalStatus = computeConstraintDerivatives();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }

  // compute secant vector
  LOCA::MultiContinuation::ExtendedVector *secantVec =
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(xVec->clone(NOX::DeepCopy));
  secantVec->update(-1.0, *prevXVec, 1.0);
  
  for (int i=0; i<numParams; i++) {
    constraints(i,0) = scaledPredictorMultiVec[i].dot(*secantVec) - 
      stepSize[i] * scaledPredictorMultiVec[i].dot(predictorMultiVec[i]);
  }

  delete secantVec;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthGroup::computeConstraintDerivatives()
{
  if (isValidConstraintDerivatives)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::NultiContinuation::ArcLengthGroup::computeConstraintDerivatives()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute predictor
  if (!isPredictor()) {
    finalStatus = computePredictor();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }

  isValidConstraintDerivatives = true;

  return finalStatus;
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraintDerivatives() const
{
  return isValidConstraintDerivatives;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::ArcLengthGroup::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::ArcLengthGroup::getConstraintDerivativesX() const
{
  return &(scaledPredictorMultiVec.getXMultiVec());
}

const NOX::Abstract::MultiVector::DenseMatrix*
LOCA::MultiContinuation::ArcLengthGroup::getConstraintDerivativesP() const
{
  return &(scaledPredictorMultiVec.getScalars());
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraintDerivativesXZero() const
{
  return false;
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraintDerivativesPZero() const
{
  return false;
}

LOCA::Extended::AbstractGroup&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
			      const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group* 
LOCA::MultiContinuation::ArcLengthGroup::clone(NOX::CopyType type) const
{
  return new LOCA::MultiContinuation::ArcLengthGroup(*this, type);
}

void
LOCA::MultiContinuation::ArcLengthGroup::notifyCompletedStep()
{
  LOCA::MultiContinuation::ExtendedGroup::notifyCompletedStep();
  isValidConstraints = false;
  isValidConstraintDerivatives = false;
}

void
LOCA::MultiContinuation::ArcLengthGroup::scalePredictor()
{
  double dpdsOld, dpdsNew;
  double thetaOld, thetaNew;
  LOCA::MultiContinuation::ExtendedVector *v, *sv;

  scaledPredictorMultiVec = predictorMultiVec;
  for (int i=0; i<numParams; i++) {
    v = 
      dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&predictorMultiVec[i]);
    sv = 
      dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&scaledPredictorMultiVec[i]);
    grpPtr->scaleVector(sv->getXVec());
    grpPtr->scaleVector(sv->getXVec());
    
    if (doArcLengthScaling) {

      // Estimate dpds
      thetaOld = theta[i];
      sv->getScalars().scale(thetaOld*thetaOld);
      dpdsOld = 1.0/sqrt(sv->dot(*v));

      if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
	cout << endl 
	     << "\t" << LOCA::Utils::fill(64, '+') << endl 
	     << "\t" << "Arc-length scaling calculation for parameter " 
	     << getContinuationParameterName(i) << ": " << endl
	     << "\t" << "Parameter component of predictor before rescaling = " 
	     << LOCA::Utils::sci(dpdsOld) << endl
	     << "\t" << "Scale factor from previous step                   = "
	     << LOCA::Utils::sci(thetaOld) << endl
	     << "\t" << "Parameter contribution to arc-length equation     = "
	     << LOCA::Utils::sci(thetaOld*dpdsOld) << endl;
      }

      // Recompute scale factor
      recalculateScaleFactor(dpdsOld, thetaOld, thetaNew);

      sv->getScalars().scale(thetaNew*thetaNew / (thetaOld*thetaOld));
      
      // Calculate new dpds using new scale factor
      dpdsNew = 1.0/sqrt(sv->dot(*v));
      
      if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
	cout << endl 
	     << "\t" << "Parameter component of predictor after rescaling  = " 
	     << LOCA::Utils::sci(dpdsNew) << endl
	     << "\t" << "New scale factor (theta)                          = "
	     << LOCA::Utils::sci(thetaNew) << endl
	     << "\t" << "Parameter contribution to arc-length equation     = "
	     << LOCA::Utils::sci(thetaNew*dpdsNew) << endl
	     << "\t" << LOCA::Utils::fill(64, '+') << endl;
      }
      
      // Rescale predictor vector
      v->scale(dpdsNew);
      sv->scale(dpdsNew);

      theta[i] = thetaNew;

      // Adjust step size scaling factor to reflect changes in 
      // arc-length parameterization
      // The first time we rescale (first continuation step) we use a 
      // different step size scale factor so that dpds*deltaS = step size 
      // provided by user
      if (isFirstRescale) {
	stepSizeScaleFactor[i] = 1.0/dpdsNew;
      }
      else
	stepSizeScaleFactor[i] = dpdsOld/dpdsNew;
    }
  }

  if (doArcLengthScaling && isFirstRescale) 
    isFirstRescale = false;
}

double
LOCA::MultiContinuation::ArcLengthGroup::computeScaledDotProduct(
			 const NOX::Abstract::Vector& x,
			 const NOX::Abstract::Vector& y) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);
  const LOCA::MultiContinuation::ExtendedVector& my = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  double val = grpPtr->computeScaledDotProduct(mx.getXVec(), my.getXVec());
  for (int i=0; i<numParams; i++)
    val += theta[i] * theta[i] * mx.getScalar(i) * my.getScalar(i);

  return val;
}

void
LOCA::MultiContinuation::ArcLengthGroup::recalculateScaleFactor(
							    double dpds,
							    double thetaOld,
							    double& thetaNew) 
{
  double g = dpds*thetaOld;

  if (g > gMax) {
    thetaNew = gGoal/dpds * sqrt( (1.0 - g*g) / (1.0 - gGoal*gGoal) ); 

    if (thetaNew < thetaMin)
      thetaNew = thetaMin;
  }
  else
    thetaNew = thetaOld;
}

void
LOCA::MultiContinuation::ArcLengthGroup::resetIsValid() {
  LOCA::MultiContinuation::ExtendedGroup::resetIsValid();
  isValidConstraints = false;
}

