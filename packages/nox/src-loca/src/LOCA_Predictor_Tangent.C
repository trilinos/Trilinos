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

#include "LOCA_Predictor_Tangent.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"

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

  // Compute derivative of residual w.r.t. parameter
  if (dfdpVecPtr == NULL)
    dfdpVecPtr = tanX.clone(NOX::ShapeCopy);
  finalStatus = underlyingGroup.computeDfDp(conParamID, *dfdpVecPtr);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // Scale dfdp by -1.0
  dfdpVecPtr->scale(-1.0);

  // Compute Jacobian
  status = underlyingGroup.computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Solve J*tanX = -df/dp
  NOX::Parameter::List& linearSolverParams = 
    LOCA::Utils::getSublist("Linear Solver");
  tanX.init(0.0);
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

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Tangent::compute(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result)
{
  string callingFunction = "LOCA::Predictor::Tangent::compute()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Get underlying group
  LOCA::MultiContinuation::AbstractGroup& underlyingGroup = 
    dynamic_cast<LOCA::MultiContinuation::AbstractGroup&>(grp.getUnderlyingGroup());

  // Get references to x, parameter components of predictor
  NOX::Abstract::MultiVector& tanX = result.getXMultiVec();
  NOX::Abstract::MultiVector::DenseMatrix& tanP = result.getScalars();

  // Get continuation parameter IDs
  const vector<int>& conParamIDs = grp.getContinuationParameterIDs();

  // Compute derivative of residual w.r.t. parameter
  NOX::Abstract::MultiVector *fdfdp = 
    xMultiVec.getXMultiVec().clone(NOX::DeepCopy);
  finalStatus = underlyingGroup.computeDfDp(conParamIDs, *fdfdp, false);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  vector<int> index_dfdp(conParamIDs.size());
  for (unsigned int i=0; i<conParamIDs.size(); i++)
    index_dfdp[i] = i+1;
  NOX::Abstract::MultiVector *dfdp = fdfdp->subView(index_dfdp);

  // Scale dfdp by -1.0
  for (unsigned int i=0; i<conParamIDs.size(); i++)
    (*dfdp)[i].scale(-1.0);

  // Compute Jacobian
  status = underlyingGroup.computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Solve J*tanX = -df/dp
  NOX::Parameter::List& linearSolverParams = 
    LOCA::Utils::getSublist("Linear Solver");
  status = underlyingGroup.applyJacobianInverseMultiVector(linearSolverParams, 
							   *dfdp, 
							   tanX);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Set parameter component equal to identity
  tanP.putScalar(0.0);
  for (unsigned int i=0; i<conParamIDs.size(); i++)
    tanP(i,i) = 1.0;

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, grp, prevXMultiVec, 
			  xMultiVec, result);

  delete fdfdp;
  delete dfdp;

  return finalStatus;
}
