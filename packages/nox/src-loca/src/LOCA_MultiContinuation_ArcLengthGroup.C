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

#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_MultiContinuation_ArcLengthConstraint.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "NOX_Parameter_List.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<NOX::Parameter::List>& continuationParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs)
  : LOCA::MultiContinuation::ExtendedGroup(global_data, topParams, 
					   continuationParams,
					   grp, pred, paramIDs),
    theta(paramIDs.size(), 1.0),
    isFirstRescale(true)
{
  Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> cons 
    = Teuchos::rcp(new LOCA::MultiContinuation::ArcLengthConstraint(
	globalData, Teuchos::rcp(this, false)));
  LOCA::MultiContinuation::ExtendedGroup::setConstraints(cons);

  double theta0 = 
    continuationParams->getParameter("Initial Scale Factor", 1.0);
  doArcLengthScaling = 
    continuationParams->getParameter("Enable Arc Length Scaling",true); 
  gGoal = 
    continuationParams->getParameter("Goal Arc Length Parameter Contribution",
				     0.5);
  gMax = 
    continuationParams->getParameter("Max Arc Length Parameter Contribution", 
				     0.0);
  thetaMin = continuationParams->getParameter("Min Scale Factor", 1.0e-3);
  
  for (int i=0; i<numParams; i++)
    theta[i] = theta0;
}

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
			 const LOCA::MultiContinuation::ArcLengthGroup& source,
			 NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type),
    theta(source.theta),
    doArcLengthScaling(source.doArcLengthScaling),
    gGoal(source.gGoal),
    gMax(source.gMax),
    thetaMin(source.thetaMin),
    isFirstRescale(source.isFirstRescale)
{
  Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ArcLengthConstraint>(conGroup->getConstraints())->setArcLengthGroup(Teuchos::rcp(this, false));
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

LOCA::Extended::MultiAbstractGroup&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
			      const LOCA::Extended::MultiAbstractGroup& source)
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

Teuchos::RefCountPtr<NOX::Abstract::Group>
LOCA::MultiContinuation::ArcLengthGroup::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ArcLengthGroup(*this, type));
}

LOCA::MultiContinuation::AbstractStrategy& 
LOCA::MultiContinuation::ArcLengthGroup::operator=(
			    const MultiContinuation::AbstractStrategy& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(source);
}

void
LOCA::MultiContinuation::ArcLengthGroup::scaleTangent()
{
  double dpdsOld, dpdsNew;
  double thetaOld, thetaNew;
  LOCA::MultiContinuation::ExtendedVector *v, *sv;

  scaledTangentMultiVec = tangentMultiVec;

  // Only scale the tangent if it is scalable
  if (predictor->isTangentScalable()) {

    for (int i=0; i<numParams; i++) {
      v = 
	dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&tangentMultiVec[i]);
      sv = 
	dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&scaledTangentMultiVec[i]);
      grpPtr->scaleVector(*(sv->getXVec()));
      grpPtr->scaleVector(*(sv->getXVec()));
      
      if (doArcLengthScaling) {
	
	// Estimate dpds
	thetaOld = theta[i];
	sv->getScalars()->scale(thetaOld*thetaOld);
	dpdsOld = 1.0/sqrt(sv->innerProduct(*v));

	if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
	  globalData->locaUtils->out() << 
	    std::endl << "\t" << 
	    globalData->locaUtils->fill(64, '+') << std::endl << "\t" << 
	    "Arc-length scaling calculation for parameter " << 
	    getContinuationParameterName(i) << ": " << std::endl << "\t" << 
	    "Parameter component of predictor before rescaling = " << 
	    globalData->locaUtils->sciformat(dpdsOld) << std::endl << "\t" << 
	    "Scale factor from previous step " << 
	    globalData->locaUtils->sciformat(thetaOld) << std::endl << "\t" << 
	    "Parameter contribution to arc-length equation     = " << 
	    globalData->locaUtils->sciformat(thetaOld*dpdsOld) << std::endl;
	}

	// Recompute scale factor
	recalculateScaleFactor(dpdsOld, thetaOld, thetaNew);
	
	sv->getScalars()->scale(thetaNew*thetaNew / (thetaOld*thetaOld));
	
	// Calculate new dpds using new scale factor
	dpdsNew = 1.0/sqrt(sv->innerProduct(*v));
	
	if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
	  globalData->locaUtils->out() << std::endl << "\t" << 
	    "Parameter component of predictor after rescaling  = " << 
	    globalData->locaUtils->sciformat(dpdsNew) << std::endl << "\t" << 
	    "New scale factor (theta)                          = " << 
	    globalData->locaUtils->sciformat(thetaNew) << std::endl << "\t" << 
	    "Parameter contribution to arc-length equation     = " << 
	    globalData->locaUtils->sciformat(thetaNew*dpdsNew) << std::endl << 
	    "\t" << globalData->locaUtils->fill(64, '+') << std::endl;
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

  double val = grpPtr->computeScaledDotProduct(*mx.getXVec(), *my.getXVec());
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
    thetaNew = gGoal/dpds * sqrt( fabs(1.0 - g*g) / fabs(1.0 - gGoal*gGoal) ); 

    if (thetaNew < thetaMin)
      thetaNew = thetaMin;
  }
  else
    thetaNew = thetaOld;
}


