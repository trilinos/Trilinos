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

#include "NOX_Parameter_List.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"
#include "LOCA_Parameter_Vector.H"

LOCA::MultiContinuation::ExtendedGroup::ExtendedGroup(
				 LOCA::MultiContinuation::AbstractGroup& g, 
				 const vector<int>& paramIDs,
				 NOX::Parameter::List& params)
  : grpPtr(&g),
    numParams(paramIDs.size()),
    xMultiVec(g.getX(), numParams+1, numParams, NOX::DeepCopy),
    fMultiVec(g.getX(), numParams+1, numParams, NOX::ShapeCopy),
    newtonMultiVec(g.getX(), numParams+1, numParams, NOX::ShapeCopy),
    gradientMultiVec(g.getX(), 1, numParams, NOX::ShapeCopy),
    predictorMultiVec(g.getX(), numParams, numParams, NOX::ShapeCopy),
    scaledPredictorMultiVec(g.getX(), numParams, numParams, NOX::ShapeCopy),
    prevXMultiVec(g.getX(), 1, numParams, NOX::ShapeCopy),
    xVec(NULL),
    fVec(NULL),
    ffMultiVec(NULL),
    dfdpMultiVec(NULL),
    newtonVec(NULL),
    gradientVec(NULL),
    predictorVec(NULL),
    scaledPredictorVec(NULL),
    prevXVec(NULL),
    borderedSolver(LOCA::Utils::getSublist("Stepper")),
    predictorManager(LOCA::Utils::getSublist("Predictor")),
    index_f(1),
    index_dfdp(numParams),
    conParamIDs(paramIDs),
    stepSize(numParams, 0.0),
    stepSizeScaleFactor(numParams, 1.0),
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false),
    isValidGradient(false),
    isValidPredictor(false),
    baseOnSecant(false)
{
  setupViews();
  for (int i=0; i<numParams; i++)
    setContinuationParameter(g.getParam(conParamIDs[i]),i);
}

LOCA::MultiContinuation::ExtendedGroup::ExtendedGroup(
				 LOCA::MultiContinuation::AbstractGroup& g, 
				 const string& paramID,
				 NOX::Parameter::List& params)
  : grpPtr(&g),
    numParams(1),
    xMultiVec(g.getX(), numParams+1, numParams, NOX::DeepCopy),
    fMultiVec(g.getX(), numParams+1, numParams, NOX::ShapeCopy),
    newtonMultiVec(g.getX(), numParams+1, numParams, NOX::ShapeCopy),
    gradientMultiVec(g.getX(), 1, numParams, NOX::ShapeCopy),
    predictorMultiVec(g.getX(), numParams, numParams, NOX::ShapeCopy),
    scaledPredictorMultiVec(g.getX(), numParams, numParams, NOX::ShapeCopy),
    prevXMultiVec(g.getX(), 1, numParams, NOX::ShapeCopy),
    xVec(NULL),
    fVec(NULL),
    ffMultiVec(NULL),
    dfdpMultiVec(NULL),
    newtonVec(NULL),
    gradientVec(NULL),
    predictorVec(NULL),
    scaledPredictorVec(NULL),
    prevXVec(NULL),
    borderedSolver(LOCA::Utils::getSublist("Stepper")),
    predictorManager(LOCA::Utils::getSublist("Predictor")),
    index_f(1),
    index_dfdp(numParams),
    conParamIDs(numParams),
    stepSize(numParams, 0.0),
    stepSizeScaleFactor(numParams, 1.0),
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false),
    isValidGradient(false),
    isValidPredictor(false),
    baseOnSecant(false)
{
  const LOCA::ParameterVector& p = g.getParams();
  conParamIDs[0] = p.getIndex(paramID);
  setupViews();
  for (int i=0; i<numParams; i++)
    setContinuationParameter(g.getParam(conParamIDs[i]),i);
}

LOCA::MultiContinuation::ExtendedGroup::ExtendedGroup(
			 const LOCA::MultiContinuation::ExtendedGroup& source,
			 NOX::CopyType type)
  : grpPtr(dynamic_cast<LOCA::MultiContinuation::AbstractGroup*>(source.grpPtr->clone())),
    numParams(source.numParams),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    gradientMultiVec(source.gradientMultiVec, type),
    predictorMultiVec(source.predictorMultiVec, type),
    scaledPredictorMultiVec(source.scaledPredictorMultiVec, type),
    prevXMultiVec(source.prevXMultiVec, type),
    xVec(NULL),
    fVec(NULL),
    ffMultiVec(NULL),
    dfdpMultiVec(NULL),
    newtonVec(NULL),
    gradientVec(NULL),
    predictorVec(NULL),
    scaledPredictorVec(NULL),
    prevXVec(NULL),
    borderedSolver(source.borderedSolver),
    predictorManager(LOCA::Utils::getSublist("Predictor")),
    index_f(1),
    index_dfdp(numParams),
    conParamIDs(source.conParamIDs),
    stepSize(source.stepSize),
    stepSizeScaleFactor(source.stepSizeScaleFactor),
    ownsGroup(true),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    isValidGradient(source.isValidGradient),
    isValidPredictor(source.isValidPredictor),
    baseOnSecant(source.baseOnSecant)
{
  setupViews();
}


LOCA::MultiContinuation::ExtendedGroup::~ExtendedGroup() 
{
  if (ownsGroup) {
    delete grpPtr;
  }
  delete ffMultiVec;
  delete dfdpMultiVec;
}

LOCA::MultiContinuation::ExtendedGroup&
LOCA::MultiContinuation::ExtendedGroup::operator=(
			 const LOCA::MultiContinuation::ExtendedGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    *grpPtr = *source.grpPtr;
    numParams = source.numParams;
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    gradientMultiVec = source.gradientMultiVec;
    predictorMultiVec = source.predictorMultiVec;
    scaledPredictorMultiVec = source.scaledPredictorMultiVec;
    prevXMultiVec = source.prevXMultiVec;
    borderedSolver = source.borderedSolver;
    conParamIDs = source.conParamIDs;
    stepSize = source.stepSize;
    stepSizeScaleFactor = source.stepSizeScaleFactor;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;
    isValidPredictor = source.isValidPredictor;
    baseOnSecant = source.baseOnSecant;

    // set up views again just to be safe
    setupViews();
  }

  return *this;
}

NOX::Abstract::Group&
LOCA::MultiContinuation::ExtendedGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(source);
}

void
LOCA::MultiContinuation::ExtendedGroup::setX(const NOX::Abstract::Vector& y)  
{
  const LOCA::MultiContinuation::ExtendedVector& my = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  grpPtr->setX( my.getXVec() );
  grpPtr->setParams(conParamIDs, my.getScalars());
  *xVec = my;

  resetIsValid();
}

void
LOCA::MultiContinuation::ExtendedGroup::computeX(
					      const NOX::Abstract::Group& g, 
					      const NOX::Abstract::Vector& d,
					      double step) 
{
  const LOCA::MultiContinuation::ExtendedGroup& mg = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(g);
  const LOCA::MultiContinuation::ExtendedVector& md = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(d);

  grpPtr->computeX(*(mg.grpPtr), md.getXVec(), step);
  xVec->update(1.0, mg.getX(), step, md, 0.0);
  grpPtr->setParams(conParamIDs, xVec->getScalars());

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  fVec->getXVec() = grpPtr->getF();
  
  // Compute constraints
  if (!isConstraints()) {
    status = computeConstraints();
  }
  fVec->getScalars().assign(getConstraints());
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  status = grpPtr->computeDfDp(conParamIDs, fMultiVec.getXMultiVec(), 
			       isValidF);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute constraint derivatives
  if (!isConstraintDerivatives()) {
    status = computeConstraintDerivatives();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  if (!isConstraintDerivativesPZero())
    dfdpMultiVec->getScalars().assign(*getConstraintDerivativesP());
  else
    dfdpMultiVec->getScalars().putScalar(0.0);

  // Set blocks in bordered solver
  borderedSolver.setIsZero(false, isConstraintDerivativesXZero(),
			   isConstraintDerivativesPZero(), false, false);
  borderedSolver.setIsContiguous(true);
  borderedSolver.setMatrixBlocks(grpPtr, &(dfdpMultiVec->getXMultiVec()), 
				 getConstraintDerivativesX(),
				 getConstraintDerivativesP());

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::computeGradient()";
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
  
  // Compute underlying gradient
  if (!grpPtr->isGradient()) {
    status = grpPtr->computeGradient();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Get grad f
  gradientVec->getXVec() = grpPtr->getGradient();

  // compute grad f + dg/dx^T * g
  if (!isConstraintDerivativesXZero())
    gradientMultiVec.update(Teuchos::TRANS, 1.0, 
		       *getConstraintDerivativesX(), 
		       getConstraints());

  // compute df/dp^T * f
  ffMultiVec->getXMultiVec().multiply(1.0, dfdpMultiVec->getXMultiVec(), 
				      gradientMultiVec.getScalars());

  // compute df/dp^T * f + dg/dp^T * g
  if (!isConstraintDerivativesPZero())
    gradientMultiVec.getScalars().multiply(Teuchos::TRANS, Teuchos::NO_TRANS,
					   1.0, *getConstraintDerivativesP(),
					   getConstraints(), 1.0);

  isValidGradient = true;

  return finalStatus;
}
   
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeNewton(
					       NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::computeNewton()";
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

  // zero out newton vec -- used as initial guess for some linear solvers
  newtonVec->init(0.0);

  status = applyJacobianInverseNewton(params);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  newtonVec->scale(-1.0);

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobian(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  NOX::Abstract::MultiVector* mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  NOX::Abstract::MultiVector* mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobian
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  // delete temporary multivectors
  delete mv_input;
  delete mv_result;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  NOX::Abstract::MultiVector* mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  NOX::Abstract::MultiVector* mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianTranspose
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianTransposeMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  // delete temporary multivectors
  delete mv_input;
  delete mv_result;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverse(
					  NOX::Parameter::List& params, 
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  NOX::Abstract::MultiVector* mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  NOX::Abstract::MultiVector* mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianInverse
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianInverseMultiVector(params, *mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  // delete temporary multivectors
  delete mv_input;
  delete mv_result;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::applyJacobianMultiVector()";
  
  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation multivectors
  const LOCA::MultiContinuation::ExtendedMultiVector& c_input = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(input);
  LOCA::MultiContinuation::ExtendedMultiVector& c_result = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::MultiVector& input_x = c_input.getXMultiVec();
  const NOX::Abstract::MultiVector::DenseMatrix& input_param = 
    c_input.getScalars();

  // Get references to x, param components of result vector
  NOX::Abstract::MultiVector& result_x = c_result.getXMultiVec();
  NOX::Abstract::MultiVector::DenseMatrix& result_param = 
    c_result.getScalars();

  // Call bordered solver apply method
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver.apply(input_x, input_param, result_x, result_param);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::applyJacobianTransposeMultiVector()";
  
  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation multivectors
  const LOCA::MultiContinuation::ExtendedMultiVector& c_input = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(input);
  LOCA::MultiContinuation::ExtendedMultiVector& c_result = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::MultiVector& input_x = c_input.getXMultiVec();
  const NOX::Abstract::MultiVector::DenseMatrix& input_param = 
    c_input.getScalars();

  // Get references to x, param components of result vector
  NOX::Abstract::MultiVector& result_x = c_result.getXMultiVec();
  NOX::Abstract::MultiVector::DenseMatrix& result_param = 
    c_result.getScalars();

  // Call bordered solver applyTranspose method
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver.applyTranspose(input_x, input_param, result_x, 
				  result_param);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverseMultiVector(
				     NOX::Parameter::List& params,
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverseMultiVector()";
  
  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation multivectors
  const LOCA::MultiContinuation::ExtendedMultiVector& c_input = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(input);
  LOCA::MultiContinuation::ExtendedMultiVector& c_result = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(result);

  // Get x, param componenets of input vector
  const NOX::Abstract::MultiVector& input_x = c_input.getXMultiVec();
  const NOX::Abstract::MultiVector::DenseMatrix& input_param = 
    c_input.getScalars();

  // Get references to x, param components of result vector
  NOX::Abstract::MultiVector& result_x = c_result.getXMultiVec();
  NOX::Abstract::MultiVector::DenseMatrix& result_param = 
    c_result.getScalars();

  // Call bordered solver applyInverse method
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver.applyInverse(params, &input_x, &input_param, 
				result_x, result_param);

  return status;
}

bool
LOCA::MultiContinuation::ExtendedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::MultiContinuation::ExtendedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::MultiContinuation::ExtendedGroup::isGradient() const 
{
  return isValidGradient;
}

bool
LOCA::MultiContinuation::ExtendedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getX() const 
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getF() const 
{
  return *fVec;
}

double
LOCA::MultiContinuation::ExtendedGroup::getNormF() const 
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getGradient() const 
{
  return *gradientVec;
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getNewton() const 
{
  return *newtonVec;
}

double
LOCA::MultiContinuation::ExtendedGroup::getNormNewtonSolveResidual() const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::MultiContinuation::ExtendedVector residual = *fVec;
  
  finalStatus = applyJacobian(*newtonVec, residual);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  residual = residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

LOCA::Extended::AbstractGroup&
LOCA::MultiContinuation::ExtendedGroup::operator=(
			      const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(source);
}

const LOCA::Continuation::AbstractGroup&
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup() const
{
  return *grpPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup()
{
  return *grpPtr;
}

int
LOCA::MultiContinuation::ExtendedGroup::getNumParams() const
{
  return numParams;
}

void
LOCA::MultiContinuation::ExtendedGroup::notifyCompletedStep()
{
  isValidPredictor = false;
  baseOnSecant = true;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiContinuation::ExtendedGroup::computePredictor()
{
  if (isValidPredictor)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType status = 
    predictorManager.compute(baseOnSecant, stepSize, *this, prevXMultiVec, 
			     xMultiVec, predictorMultiVec);

  scalePredictor();

  isValidPredictor = true;
  return status;
}

const LOCA::MultiContinuation::ExtendedVector&
LOCA::MultiContinuation::ExtendedGroup::getPredictorDirection(int i) const 
{
  return dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(predictorMultiVec[i]);
}

const LOCA::MultiContinuation::ExtendedMultiVector&
LOCA::MultiContinuation::ExtendedGroup::getPredictorDirections() const 
{
  return predictorMultiVec;
}

void
LOCA::MultiContinuation::ExtendedGroup::setPredictorDirection(
			    const LOCA::MultiContinuation::ExtendedVector& v, 
			    int i)
{
  predictorMultiVec[i] = v;
  isValidPredictor = true;
}

void
LOCA::MultiContinuation::ExtendedGroup::setPredictorDirections(
			const LOCA::MultiContinuation::ExtendedMultiVector& v)
{
  predictorMultiVec = v;
  isValidPredictor = true;
}

bool
LOCA::MultiContinuation::ExtendedGroup::isPredictor() const
{
  return isValidPredictor;
}

void 
LOCA::MultiContinuation::ExtendedGroup::resetPredictor(
				       NOX::Parameter::List& predictorParams)
{
  baseOnSecant = false;
  predictorManager.reset(predictorParams);
}

void
LOCA::MultiContinuation::ExtendedGroup::scalePredictor()
{
  LOCA::MultiContinuation::ExtendedVector *v;

  scaledPredictorMultiVec = predictorMultiVec;
  for (int i=0; i<numParams; i++) {
    v = 
      dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&scaledPredictorMultiVec[i]);
    grpPtr->scaleVector(v->getXVec());
    grpPtr->scaleVector(v->getXVec());
  }
}

void
LOCA::MultiContinuation::ExtendedGroup::setPrevX(
					     const NOX::Abstract::Vector& y) 
{
  *prevXVec = y;
}

const LOCA::MultiContinuation::ExtendedVector&
LOCA::MultiContinuation::ExtendedGroup::getPrevX() const 
{
  return *prevXVec;
}

void
LOCA::MultiContinuation::ExtendedGroup::setStepSize(double deltaS, int i) 
{
  stepSize[i] = deltaS;
}

double
LOCA::MultiContinuation::ExtendedGroup::getStepSize(int i) const 
{
  return stepSize[i];
}

void
LOCA::MultiContinuation::ExtendedGroup::setContinuationParameter(double val,
								 int i) 
{
  grpPtr->setParam(conParamIDs[i], val);
  xVec->getScalar(i) = val;
}

double
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameter(int i) const 
{
  return grpPtr->getParam(conParamIDs[i]);
}

int
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterID(int i) const
{
  return conParamIDs[i];
}

const vector<int>&
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterIDs() const
{
  return conParamIDs;
}

string
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterName(
								  int i) const
{
  const LOCA::ParameterVector& p = grpPtr->getParams();
  return p.getLabel(conParamIDs[i]);
}

double
LOCA::MultiContinuation::ExtendedGroup::getStepSizeScaleFactor(int i) const 
{
  return stepSizeScaleFactor[i];
}

void
LOCA::MultiContinuation::ExtendedGroup::printSolution() const
{
  for (int i=0; i<numParams; i++)
    grpPtr->printSolution(getContinuationParameter(i));
}

double
LOCA::MultiContinuation::ExtendedGroup::computeScaledDotProduct(
			 const NOX::Abstract::Vector& x,
			 const NOX::Abstract::Vector& y) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);
  const LOCA::MultiContinuation::ExtendedVector& my = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  double val = grpPtr->computeScaledDotProduct(mx.getXVec(), my.getXVec());
  for (int i=0; i<numParams; i++)
    val += mx.getScalar(i) * my.getScalar(i);

  return val;
}

int
LOCA::MultiContinuation::ExtendedGroup::projectToDrawDimension() const
{
  return numParams + grpPtr->projectToDrawDimension();
}

void
LOCA::MultiContinuation::ExtendedGroup::projectToDraw(
			    const LOCA::MultiContinuation::ExtendedVector& x, 
			    double *px) const
{
  // first numParams components are the parameters
  for (int i=0; i<numParams; i++)
    px[i] = x.getScalar(i);

  // fill remaining solution components
  grpPtr->projectToDraw(x.getXVec(), px+numParams);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverseNewton(
						NOX::Parameter::List& params) 
{
  // This method is specialized to the Newton solve where the right-hand-side
  // is f.  We take advantage of the fact that f and df/dp are in a 
  // contiguous multivector

  string callingFunction = 
    "LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverseNewton()";
  
  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Get x, param components of f vector (we only want the parameter 
  // components of f, not df/dp)
  const NOX::Abstract::MultiVector& f_x = fMultiVec.getXMultiVec();
  const NOX::Abstract::MultiVector::DenseMatrix& f_p = 
    ffMultiVec->getScalars();

  // Get references to x, param components of newton vector
  NOX::Abstract::MultiVector& newton_x = newtonMultiVec.getXMultiVec();
  NOX::Abstract::MultiVector::DenseMatrix& newton_p = 
    newtonVec->getScalars();

  // Call bordered solver applyInverse method
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver.applyInverse(params, &f_x, &f_p, newton_x, newton_p);

  return status;
}

void
LOCA::MultiContinuation::ExtendedGroup::resetIsValid() {
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
}

void
LOCA::MultiContinuation::ExtendedGroup::setupViews()
{
  index_f[0] = 0;
  for (int i=0; i<numParams; i++)
    index_dfdp[i] = i+1;

  
  xVec = dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&xMultiVec[0]);
  fVec = dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&fMultiVec[0]);
  newtonVec = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&newtonMultiVec[0]);
  gradientVec = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&gradientMultiVec[0]);
  predictorVec = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&predictorMultiVec[0]);
  scaledPredictorVec = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&scaledPredictorMultiVec[0]);
  prevXVec = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector*>(&prevXMultiVec[0]);

  if (ffMultiVec != NULL)
    delete ffMultiVec;
  ffMultiVec = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector*>(fMultiVec.subView(index_f));

  if (dfdpMultiVec != NULL)
    delete dfdpMultiVec;
  dfdpMultiVec = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector*>(fMultiVec.subView(index_dfdp));

  

}
