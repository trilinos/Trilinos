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

#include "LOCA_Predictor_Generic.H"
#include "LOCA_Continuation_ExtendedVector.H"
#include "LOCA_Continuation_ExtendedGroup.H"

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
