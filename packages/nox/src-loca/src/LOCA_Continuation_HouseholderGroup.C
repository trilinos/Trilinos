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

#include "LOCA_Continuation_HouseholderGroup.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"

LOCA::Continuation::HouseholderGroup::HouseholderGroup(
				 LOCA::Continuation::AbstractGroup& g,
				 int paramID,
				 NOX::Parameter::List& params)
  : LOCA::Continuation::ArcLengthGroup(g, paramID, params), 
    houseVec(g.getX(), 0.0),
    beta(0.0)
{
}

LOCA::Continuation::HouseholderGroup::HouseholderGroup(
				  LOCA::Continuation::AbstractGroup& g,
				  string paramID,
				  NOX::Parameter::List& params)
  : LOCA::Continuation::ArcLengthGroup(g, paramID, params), 
    houseVec(g.getX(), 0.0),
    beta(0.0)
{
}

LOCA::Continuation::HouseholderGroup::HouseholderGroup(
		   const LOCA::Continuation::HouseholderGroup& source,
		   NOX::CopyType type)
  :  LOCA::Continuation::ArcLengthGroup(source, type), 
     houseVec(source.houseVec, type),
     beta(source.beta)
{
}


LOCA::Continuation::HouseholderGroup::~HouseholderGroup() 
{
}

LOCA::Continuation::ArcLengthGroup&
LOCA::Continuation::HouseholderGroup::operator=(
			      const LOCA::Continuation::ArcLengthGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Continuation::HouseholderGroup&>(source);
}

LOCA::Continuation::ExtendedGroup&
LOCA::Continuation::HouseholderGroup::operator=(
			      const LOCA::Continuation::ExtendedGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Continuation::HouseholderGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Continuation::HouseholderGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::Continuation::HouseholderGroup&>(source);
}

LOCA::Continuation::HouseholderGroup&
LOCA::Continuation::HouseholderGroup::operator=(
		   const LOCA::Continuation::HouseholderGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    LOCA::Continuation::ArcLengthGroup::operator=(source);

    // Copy values
    houseVec = source.houseVec;
    beta = source.beta;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Continuation::HouseholderGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Continuation::HouseholderGroup(*this);
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::HouseholderGroup::computeF() 
{
  // Call computeF() of arclength group
  NOX::Abstract::Group::ReturnType finalStatus = 
    LOCA::Continuation::ArcLengthGroup::computeF();
  
  // Set parameter component to zero
  fVec.getParam() = 0.0;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::HouseholderGroup::computeNewton(
					       NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Continuation::HouseholderGroup::computeNewton()";
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

  // Compute Householder vector
  computeHouseholderVector();

  // zero out newton vec -- used as initial guess for some linear solvers
  newtonVec.init(0.0);

  // Get x component of fVec, newtonVec, houseVec
  const NOX::Abstract::Vector& fVec_x = fVec.getXVec();
  const NOX::Abstract::Vector& houseVec_x = houseVec.getXVec();
  NOX::Abstract::Vector& newtonVec_x = newtonVec.getXVec();
  
  // Get param component of houseVec, newtonVec
  const double& houseVec_p = houseVec.getParam();
  double& newtonVec_p = newtonVec.getParam();

  // Solve continuation equations using Householder projection
  status = grpPtr->applyHouseholderJacobianInverse(params, fVec_x, 
						   *derivResidualParamPtr,
						   houseVec_x, houseVec_p, 
						   beta, newtonVec_x,
						   newtonVec_p);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Change sign of newton vector
  newtonVec.scale(-1.0);

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    double r = computeScaledDotProduct(newtonVec, predictorVec);
    cout << "\n\tScaled component of Newton vector in direction of "
	 << "predictor:  " << r << endl;
  }

  isValidNewton = true;

  return finalStatus;
}

void
LOCA::Continuation::HouseholderGroup::computeHouseholderVector()
{

  // Copy predictor in Householder vector and scale (twice)
  houseVec = predictorVec;
  scaleVector(houseVec);
  scaleVector(houseVec);

  // Get x, param components of Householder vector
  NOX::Abstract::Vector& houseX = houseVec.getXVec();
  double& houseP = houseVec.getParam();

  double sigma = houseX.dot(houseX);

  if (sigma == 0.0) {
    beta = 0.0;
    houseP = 1.0;
  }
  else {
    double mu = sqrt(houseP*houseP + sigma);
    if (houseP <= 0.0)
      houseP = houseP - mu;
    else
      houseP = -sigma / (houseP + mu);
    beta = 2.0*houseP*houseP/(sigma + houseP*houseP);
    houseVec.scale(1.0/houseP);
  }

  return;
}

// void
// LOCA::Continuation::HouseholderGroup::computeHouseholderVector()
// {

//   // Scale predictor
//   LOCA::Continuation::ExtendedVector *scaledPredictor = 
//     dynamic_cast<LOCA::Continuation::ExtendedVector*>(predictorVec.clone(NOX::DeepCopy));
//   scaleVector(*scaledPredictor);
//   scaleVector(*scaledPredictor);
  
//   // Get x, param components of predictor vector
//   const NOX::Abstract::Vector& tanX = scaledPredictor->getXVec();
//   const double& tanP = scaledPredictor->getParam();

//   // Get x, param components of Householder vector
//   NOX::Abstract::Vector& houseX = houseVec.getXVec();
//   double& houseP = houseVec.getParam();

//   double sigma = tanX.dot(tanX);
//   houseX = tanX;

//   if (sigma == 0.0) {
//     beta = 0.0;
//     houseP = 1.0;
//   }
//   else {
//     double mu = sqrt(tanP*tanP + sigma);
//     if (tanP <= 0.0)
//       houseP = tanP - mu;
//     else
//       houseP = -sigma / (tanP + mu);
//     beta = 2.0*houseP*houseP/(sigma + houseP*houseP);
//     houseVec.scale(1.0/houseP);
//   }

//   delete scaledPredictor;

//   return;
// }

void
LOCA::Continuation::HouseholderGroup::scaleVector(
				 LOCA::Continuation::ExtendedVector& x) const
{
  grpPtr->scaleVector(x.getXVec());
  x.getParam() *= theta;
}



