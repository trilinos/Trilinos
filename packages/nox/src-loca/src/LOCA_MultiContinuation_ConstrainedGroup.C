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

#include "NOX_Parameter_List.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSystem_AbstractStrategy.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"
#include "LOCA_Parameter_Vector.H"

LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup(
       const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<NOX::Parameter::List>& conParams,
       const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& g,
       const Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>& constraints,
       const vector<int>& paramIDs)
  : globalData(global_data),
    parsedParams(topParams),
    constraintParams(conParams),
    grpPtr(g),
    constraintsPtr(constraints),
    numParams(paramIDs.size()),
    xMultiVec(globalData, g->getX(), numParams+1, numParams, NOX::DeepCopy),
    fMultiVec(globalData, g->getX(), numParams+1, numParams, NOX::ShapeCopy),
    newtonMultiVec(globalData, g->getX(), numParams+1, numParams, 
		   NOX::ShapeCopy),
    gradientMultiVec(globalData, g->getX(), 1, numParams, NOX::ShapeCopy),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    gradientVec(),
    borderedSolver(),
    index_f(1),
    index_dfdp(numParams),
    constraintParamIDs(paramIDs),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false),
    isValidGradient(false)
{
  // Set up multi-vector views
  setupViews(); 

  // Set parameters in solution vector
  for (int i=0; i<numParams; i++)
    xVec->getScalar(i) = grpPtr->getParam(constraintParamIDs[i]);

  // Set parameters and solution vector in constraints
  constraintsPtr->setParams(constraintParamIDs, *xVec->getScalars());
  constraintsPtr->setX(*(xVec->getXVec()));

  // Instantiate bordered solver
  borderedSolver = globalData->locaFactory->createBorderedSystemStrategy(
				   parsedParams,
				   constraintParams);
}

LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup(
		     const LOCA::MultiContinuation::ConstrainedGroup& source,
		     NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    constraintParams(source.constraintParams),
    grpPtr(Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::AbstractGroup>(source.grpPtr->clone(type))),
    constraintsPtr(source.constraintsPtr->clone(type)),
    numParams(source.numParams),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    gradientMultiVec(source.gradientMultiVec, type),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    gradientVec(),
    borderedSolver(source.borderedSolver),
    index_f(1),
    index_dfdp(numParams),
    constraintParamIDs(source.constraintParamIDs),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    isValidGradient(source.isValidGradient)
{
  // Set up multi-vector views
  setupViews();

  // Instantiate bordered solver
  borderedSolver = globalData->locaFactory->createBorderedSystemStrategy(
				   parsedParams,
				   constraintParams);

  if (type == NOX::ShapeCopy) {
    isValidF = false;
    isValidJacobian = false;
    isValidNewton = false;
    isValidGradient = false;
  }
}


LOCA::MultiContinuation::ConstrainedGroup::~ConstrainedGroup() 
{
}

void
LOCA::MultiContinuation::ConstrainedGroup::setConstraintParameter(int i,
								  double val) 
{
  grpPtr->setParam(constraintParamIDs[i],val);
  xVec->getScalar(i) = val;
  constraintsPtr->setParam(constraintParamIDs[i],i);

  resetIsValid();
}

double
LOCA::MultiContinuation::ConstrainedGroup::getConstraintParameter(int i) const
{
  return grpPtr->getParam(constraintParamIDs[i]);
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ConstrainedGroup::getGroup() const
{
  return grpPtr;
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::ConstrainedGroup::getConstraints() const
{
  return constraintsPtr;
}

const vector<int>&
LOCA::MultiContinuation::ConstrainedGroup::getConstraintParamIDs() const
{
  return constraintParamIDs;
}

NOX::Abstract::Group&
LOCA::MultiContinuation::ConstrainedGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RefCountPtr<NOX::Abstract::Group>
LOCA::MultiContinuation::ConstrainedGroup::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new ConstrainedGroup(*this, type));
}

void
LOCA::MultiContinuation::ConstrainedGroup::setX(
					     const NOX::Abstract::Vector& y)  
{
  const LOCA::MultiContinuation::ExtendedVector& my = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  grpPtr->setX( *(my.getXVec()) );
  grpPtr->setParamsMulti(constraintParamIDs, *my.getScalars());
  *xVec = my;
  constraintsPtr->setX( *(my.getXVec()) );
  constraintsPtr->setParams(constraintParamIDs, *my.getScalars());

  resetIsValid();
}

void
LOCA::MultiContinuation::ConstrainedGroup::computeX(
					      const NOX::Abstract::Group& g, 
					      const NOX::Abstract::Vector& d,
					      double step) 
{
  const LOCA::MultiContinuation::ConstrainedGroup& mg = 
    dynamic_cast<const LOCA::MultiContinuation::ConstrainedGroup&>(g);
  const LOCA::MultiContinuation::ExtendedVector& md = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(d);

  grpPtr->computeX(*(mg.grpPtr), *(md.getXVec()), step);
  xVec->update(1.0, mg.getX(), step, md, 0.0);
  grpPtr->setParamsMulti(constraintParamIDs, *xVec->getScalars());
  constraintsPtr->setX( *(xVec->getXVec()) );
  constraintsPtr->setParams(constraintParamIDs, *xVec->getScalars());

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::computeF()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }
  *(fVec->getXVec()) = grpPtr->getF();
  
  // Compute constraints
  if (!constraintsPtr->isConstraints()) {
    status = constraintsPtr->computeConstraints();
  }
  fVec->getScalars()->assign(constraintsPtr->getConstraints());
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  status = grpPtr->computeDfDpMulti(constraintParamIDs, 
				    *fMultiVec.getXMultiVec(), 
				    isValidF);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							    finalStatus,
							    callingFunction);

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }

  // Compute constraint derivatives
  if (!constraintsPtr->isDX()) {
    status = constraintsPtr->computeDX();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }
  status = 
    constraintsPtr->computeDP(constraintParamIDs,
			      *fMultiVec.getScalars(),
			      isValidF);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							    finalStatus,
							    callingFunction);

  // Set blocks in bordered solver
  borderedSolver->setMatrixBlocks(
			 grpPtr, 
			 dfdpMultiVec->getXMultiVec(), 
			 constraintsPtr,
			 dfdpMultiVec->getScalars());

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::computeGradient()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }
  
  // Make sure Jacobian is valid
  if (!isJacobian()) {
    status = computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }
  
  // Compute underlying gradient
  if (!grpPtr->isGradient()) {
    status = grpPtr->computeGradient();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }

  // Get grad f
  *gradientVec->getXVec() = grpPtr->getGradient();

  // compute grad f + dg/dx^T * g
  constraintsPtr->addDX(Teuchos::TRANS, 1.0, 
			constraintsPtr->getConstraints(),
			1.0, 
			*gradientMultiVec.getXMultiVec());

  // compute df/dp^T * f
  ffMultiVec->getXMultiVec()->multiply(1.0, *dfdpMultiVec->getXMultiVec(), 
				       *gradientMultiVec.getScalars());

  // compute df/dp^T * f + dg/dp^T * g
  gradientMultiVec.getScalars()->multiply(
			       Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, 
			       *dfdpMultiVec->getScalars(),
			       constraintsPtr->getConstraints(), 1.0);

  isValidGradient = true;

  return finalStatus;
}
   
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeNewton(
					       NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::computeNewton()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }
  
  // Make sure Jacobian is valid
  if (!isJacobian()) {
    status = computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							      finalStatus,
							      callingFunction);
  }

  // zero out newton vec -- used as initial guess for some linear solvers
  newtonVec->init(0.0);

  status = applyJacobianInverseNewton(params);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							    finalStatus,
							    callingFunction);

  newtonVec->scale(-1.0);

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobian(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobian
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianTranspose
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianTransposeMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverse(
					  NOX::Parameter::List& params, 
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianInverse
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianInverseMultiVector(params, *mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::applyJacobianMultiVector()";
  
  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
					    "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation multivectors
  const LOCA::MultiContinuation::ExtendedMultiVector& c_input = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(input);
  LOCA::MultiContinuation::ExtendedMultiVector& c_result = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(result);

  // Get x, param componenets of input vector
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> input_x = 
    c_input.getXMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_x = 
    c_result.getXMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    c_result.getScalars();

  // Call bordered solver apply method
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver->apply(*input_x, *input_param, *result_x, *result_param);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeMultiVector()";
  
  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
					    "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation multivectors
  const LOCA::MultiContinuation::ExtendedMultiVector& c_input = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(input);
  LOCA::MultiContinuation::ExtendedMultiVector& c_result = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(result);

  // Get x, param componenets of input vector
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> input_x = 
    c_input.getXMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_x = 
    c_result.getXMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    c_result.getScalars();

  // Call bordered solver applyTranspose method
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver->applyTranspose(*input_x, *input_param, *result_x, 
				   *result_param);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseMultiVector(
				     NOX::Parameter::List& params,
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseMultiVector()";
  
  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
					    "Called with invalid Jacobian!");
  }

  // Cast inputs to continuation multivectors
  const LOCA::MultiContinuation::ExtendedMultiVector& c_input = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(input);
  LOCA::MultiContinuation::ExtendedMultiVector& c_result = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(result);

  // Get x, param componenets of input vector
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> input_x = 
    c_input.getXMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_x = 
    c_result.getXMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    c_result.getScalars();

  // Call bordered solver applyInverse method
  borderedSolver->setIsContiguous(false);
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver->applyInverse(params, input_x.get(), input_param.get(), 
				 *result_x, *result_param);

  return status;
}

bool
LOCA::MultiContinuation::ConstrainedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::MultiContinuation::ConstrainedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::MultiContinuation::ConstrainedGroup::isGradient() const 
{
  return isValidGradient;
}

bool
LOCA::MultiContinuation::ConstrainedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::MultiContinuation::ConstrainedGroup::getX() const 
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ConstrainedGroup::getF() const 
{
  return *fVec;
}

double
LOCA::MultiContinuation::ConstrainedGroup::getNormF() const 
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ConstrainedGroup::getGradient() const 
{
  return *gradientVec;
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ConstrainedGroup::getNewton() const 
{
  return *newtonVec;
}

double
LOCA::MultiContinuation::ConstrainedGroup::getNormNewtonSolveResidual() const 
{
  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::MultiContinuation::ExtendedVector residual = *fVec;
  
  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual = residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup()
{
  return grpPtr;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseNewton(
						NOX::Parameter::List& params) 
{
  // This method is specialized to the Newton solve where the right-hand-side
  // is f.  We take advantage of the fact that f and df/dp are in a 
  // contiguous multivector

  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseNewton()";
  
  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
					    "Called with invalid Jacobian!");
  }

  // Get x, param components of f vector (we only want the parameter 
  // components of f, not df/dp)
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> f_x = 
    fMultiVec.getXMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> f_p = 
    ffMultiVec->getScalars();

  // Get references to x, param components of newton vector
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> newton_x = 
    newtonMultiVec.getXMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> newton_p = 
    newtonVec->getScalars();

  // Call bordered solver applyInverse method
  borderedSolver->setIsContiguous(true);
  NOX::Abstract::Group::ReturnType status = 
    borderedSolver->applyInverse(params, f_x.get(), f_p.get(), *newton_x, 
				 *newton_p);

  return status;
}

void
LOCA::MultiContinuation::ConstrainedGroup::copy(
					      const NOX::Abstract::Group& src) 
{

  const LOCA::MultiContinuation::ConstrainedGroup& source = 
    dynamic_cast<const LOCA::MultiContinuation::ConstrainedGroup&>(src);

  // Protect against A = A
  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    constraintParams = source.constraintParams;
    grpPtr->copy(*source.grpPtr);
    constraintsPtr->copy(*source.constraintsPtr);
    numParams = source.numParams;
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    gradientMultiVec = source.gradientMultiVec;
    index_f = source.index_f;
    index_dfdp = source.index_dfdp;
    constraintParamIDs = source.constraintParamIDs;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;

    // set up views again just to be safe
    setupViews();

    // Instantiate bordered solver
    borderedSolver = 
      globalData->locaFactory->createBorderedSystemStrategy(
				   parsedParams,
				   constraintParams);
  }
}

void
LOCA::MultiContinuation::ConstrainedGroup::setParamsMulti(
		     const vector<int>& paramIDs, 
		     const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  constraintsPtr->setParams(paramIDs, vals);

  for (unsigned int i=0; i<paramIDs.size(); i++)
    for (unsigned int j=0; j<constraintParamIDs.size(); j++)
      if (paramIDs[i] == constraintParamIDs[j])
	xVec->getScalar(j) = vals(i,0);

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeDfDpMulti(
					     const vector<int>& paramIDs, 
					     NOX::Abstract::MultiVector& dfdp, 
					     bool isValid_F)
{
  string callingFunction = 
    "LOCA::MultiContinuation::ConstrainedGroup::computeDfDpMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast result to constraint vector
  LOCA::MultiContinuation::ExtendedMultiVector& c_dfdp = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(dfdp);

  // Compute df/dp
  status = grpPtr->computeDfDpMulti(paramIDs, *c_dfdp.getXMultiVec(), 
				    isValid_F);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							    finalStatus,
							    callingFunction);

  // Compute dg/dp
  status = constraintsPtr->computeDP(paramIDs, 
				     *c_dfdp.getScalars(), 
				     isValid_F);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							    finalStatus,
							    callingFunction);

  return finalStatus;
}

void
LOCA::MultiContinuation::ConstrainedGroup::projectToDraw(
					      const NOX::Abstract::Vector& x,
					      double *px) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);

  grpPtr->projectToDraw(*mx.getXVec(), px);
  for (int i=0; i<numParams; i++) {
    px[grpPtr->projectToDrawDimension()+i] = mx.getScalar(i);
  }
}

int
LOCA::MultiContinuation::ConstrainedGroup::projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension() + numParams;
}

void
LOCA::MultiContinuation::ConstrainedGroup::setParams(
					       const LOCA::ParameterVector& p)
{
  grpPtr->setParams(p);
  for (int i=0; i<p.length(); i++)
    constraintsPtr->setParam(i, p[i]);
  for (int i=0; i<numParams; i++)
    xVec->getScalar(i) = p[constraintParamIDs[i]];

  resetIsValid();
}

void
LOCA::MultiContinuation::ConstrainedGroup::setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
  constraintsPtr->setParam(paramID, val);

  for (unsigned int j=0; j<constraintParamIDs.size(); j++)
    if (paramID == constraintParamIDs[j])
      xVec->getScalar(j) = val;

  resetIsValid();
}

void
LOCA::MultiContinuation::ConstrainedGroup::setParam(string paramID, double val)
{
  const LOCA::ParameterVector& p = grpPtr->getParams();
  int id = p.getIndex(paramID);
  setParam(id, val);
}

const LOCA::ParameterVector&
LOCA::MultiContinuation::ConstrainedGroup::getParams() const
{
  return grpPtr->getParams();
}

double
LOCA::MultiContinuation::ConstrainedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::MultiContinuation::ConstrainedGroup::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::MultiContinuation::ConstrainedGroup::computeScaledDotProduct(
				       const NOX::Abstract::Vector& a,
				       const NOX::Abstract::Vector& b) const
{
  const LOCA::MultiContinuation::ExtendedVector& ma = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(a);
  const LOCA::MultiContinuation::ExtendedVector& mb = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(b);

  double val = grpPtr->computeScaledDotProduct(*ma.getXVec(), *mb.getXVec());
  for (int i=0; i<numParams; i++) {
    val += ma.getScalar(i) * mb.getScalar(i);
  }

  return val;
}

void
LOCA::MultiContinuation::ConstrainedGroup::printSolution(
						 const double conParam) const
{
  printSolution(*xVec, conParam);
}

void
LOCA::MultiContinuation::ConstrainedGroup::printSolution(
					      const NOX::Abstract::Vector& x,
					      const double conParam) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "LOCA::MultiContinuation::ConstrainedGroup::printSolution\n";

    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for conParam = " << 
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(*mx.getXVec(), conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << "\tPrinting constraint parameters\n";
    mx.getScalars()->print(globalData->locaUtils->out());
  }
}

void
LOCA::MultiContinuation::ConstrainedGroup::scaleVector(
					       NOX::Abstract::Vector& x) const
{
  LOCA::MultiContinuation::ExtendedVector& mx = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector&>(x);

  grpPtr->scaleVector(*mx.getXVec());
}

void
LOCA::MultiContinuation::ConstrainedGroup::resetIsValid() {
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
}

void
LOCA::MultiContinuation::ConstrainedGroup::setupViews()
{
  index_f[0] = 0;
  for (int i=0; i<numParams; i++)
    index_dfdp[i] = i+1;
  
  xVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(xMultiVec.getVector(0),true);
  fVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(fMultiVec.getVector(0),true);
  newtonVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(newtonMultiVec.getVector(0),true);
  gradientVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(gradientMultiVec.getVector(0),true);

  ffMultiVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(fMultiVec.subView(index_f),true);

  dfdpMultiVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(fMultiVec.subView(index_dfdp),true);

}
