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

#include "LOCA_Predictor_Tangent.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::Predictor::Tangent::Tangent(NOX::Parameter::List& params) :
  dfdpVecPtr(NULL)
{
}

LOCA::Predictor::Tangent::~Tangent()
{
  if (dfdpVecPtr != NULL)
    delete dfdpVecPtr;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Tangent::reset(NOX::Parameter::List& params) 
{
  if (dfdpVecPtr != NULL) {
    delete dfdpVecPtr;
    dfdpVecPtr = NULL;
  }
  return LOCA::Predictor::Generic::reset(params);
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Tangent::compute(bool baseOnSecant, double stepSize,
				  LOCA::Continuation::ExtendedGroup& prevGroup,
				  LOCA::Continuation::ExtendedGroup& curGroup,
				  LOCA::Continuation::ExtendedVector& result) 
{
  string callingFunction = "LOCA::Predictor::Tangent::compute()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Get references to x, parameter components of predictor
  NOX::Abstract::Vector& tanX = result.getXVec();
  double& tanP = result.getParam();

  // Get underlying group
  LOCA::Continuation::AbstractGroup& underlyingGroup = 
    curGroup.getUnderlyingGroup();

  // Get continuation parameter ID
  int conParamID = curGroup.getContinuationParameterID();

  // Compute Jacobian
  finalStatus = underlyingGroup.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // Compute derivative of residual w.r.t. parameter
  if (dfdpVecPtr == NULL)
    dfdpVecPtr = tanX.clone(NOX::ShapeCopy);
  status = underlyingGroup.computeDfDp(conParamID, *dfdpVecPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Scale dfdp by -1.0
  dfdpVecPtr->scale(-1.0);
  
  // Solve J*tanX = -df/dp
  NOX::Parameter::List& linearSolverParams = 
    LOCA::Utils::getSublist("Linear Solver");
  status = underlyingGroup.applyJacobianInverse(linearSolverParams, 
						*dfdpVecPtr, 
						tanX);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Set parameter component equal to 1
  tanP = 1.0;

  // Rescale predictor
  curGroup.scalePredictor(result);

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, prevGroup, curGroup, result);
  
  // Set predictor in continuation group
  curGroup.setPredictorDirection(result);

  return finalStatus;
}
