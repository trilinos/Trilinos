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

#include "LOCA_Predictor_Generic.H"
#include "LOCA_Continuation_ExtendedVector.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"

LOCA::Predictor::Generic::Generic()
  : secantVecPtr(NULL) {}

LOCA::Predictor::Generic::~Generic()
{
  if (secantVecPtr != NULL)
    delete secantVecPtr;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Generic::reset(NOX::Parameter::List& params)
{
  if (secantVecPtr != NULL) {
    delete secantVecPtr;
    secantVecPtr = NULL;
  }

  return NOX::Abstract::Group::Ok;
}

void
LOCA::Predictor::Generic::setPredictorOrientation(
				bool baseOnSecant, double ds,
				LOCA::Continuation::ExtendedGroup& prevGroup,
				LOCA::Continuation::ExtendedGroup& curGroup,
				LOCA::Continuation::ExtendedVector& result) 
{
  // If orientation is not based on a secant vector (i.e., first or last
  // steps in a continuation run) make parameter component of predictor
  // positive
  if (!baseOnSecant) {
    if (result.getParam() < 0.0)
      result.scale(-1.0);
    return;
  }

  // Compute secant vector
  if (secantVecPtr == NULL) {
    secantVecPtr = 
      dynamic_cast<LOCA::Continuation::ExtendedVector*>(curGroup.getX().clone(NOX::DeepCopy));
  }
  else {
    *secantVecPtr = curGroup.getX();
  }
  secantVecPtr->update(-1.0, prevGroup.getX(), 1.0);

  if (curGroup.computeScaledDotProduct(*secantVecPtr, result)*ds < 0.0)
    result.scale(-1.0);
}

void
LOCA::Predictor::Generic::setPredictorOrientation(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result)
{
  // If orientation is not based on a secant vector (i.e., first or last
  // steps in a continuation run) make parameter component of predictor
  // positive
  if (!baseOnSecant) {
    for (int i=0; i<result.numVectors(); i++) 
      if (result.getScalar(i,i) < 0.0)
	result[i].scale(-1.0);
    return;
  }

  LOCA::MultiContinuation::ExtendedVector* secantVecPtr = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(xMultiVec[0].clone(NOX::DeepCopy));
  secantVecPtr->update(-1.0, prevXMultiVec[0], 1.0);

  for (int i=0; i<result.numVectors(); i++)
    if (secantVecPtr->dot(result[i])*stepSize[i] < 0.0)
      result[i].scale(-1.0);

  delete secantVecPtr;
}
