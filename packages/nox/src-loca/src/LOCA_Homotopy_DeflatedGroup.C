// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_Homotopy_DeflatedGroup.H"

#include "Teuchos_ParameterList.hpp"
#include "LOCA_Homotopy_AbstractGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"
#include "LOCA_BorderedSolver_JacobianOperator.H"

LOCA::Homotopy::DeflatedGroup::
DeflatedGroup(
       const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const Teuchos::RCP<Teuchos::ParameterList>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& hParams,
       const Teuchos::RCP<LOCA::Homotopy::AbstractGroup>& g,
       const Teuchos::RCP<const NOX::Abstract::Vector>& start_vec,
       const std::vector< Teuchos::RCP<const NOX::Abstract::Vector> >& prev_solns,
       const double identity_sign)
  : globalData(global_data),
    parsedParams(),
    homotopyParams(hParams),
    grpPtr(g),
    bordered_grp(),
    xMultiVec(globalData, g->getX(), 1, 1, NOX::DeepCopy),
    fMultiVec(globalData, g->getX(), 1, 1, NOX::ShapeCopy),
    newtonMultiVec(globalData, g->getX(), 1, 1, NOX::ShapeCopy),
    gradientMultiVec(globalData, g->getX(), 1, 1, NOX::ShapeCopy),
    xVec(),
    fVec(),
    newtonVec(),
    gradientVec(),
    startVec(start_vec),
    identitySign(identity_sign),
    solns(prev_solns),
    distVec(startVec->clone(NOX::ShapeCopy)),
    totalDistMultiVec(startVec->createMultiVector(1, NOX::ShapeCopy)),
    underlyingF(startVec->createMultiVector(1, NOX::ShapeCopy)),
    jacOp(),
    borderedSolver(),
    minusOne(),
    numSolns(solns.size()),
    distances(numSolns),
    distProd(0.0),
    index_f(1),
    paramVec(grpPtr->getParams()),
    conParam(0.0),
    conParamID(-1),
    conParamLabel("Homotopy Continuation Parameter"),
    augmentJacForHomotopyNotImplemented(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false),
    isValidGradient(false),
    isBordered(false)
{
  // Set up multi-vector views
  setupViews(); 

  // Set the homotopy parameter as a parameter in the ParameterVector.
  // This will allow us to do an invasive homotopy since the app can now 
  // get access to the homotopy continuation parameter.
  paramVec.addParameter(conParamLabel, conParam);
  grpPtr->setParams(paramVec);

  // Set up paramID
  conParamID = paramVec.getIndex(conParamLabel);

  setStepperParameters(*topParams);

  // Parse parameter list
  parsedParams = Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
  parsedParams->parseSublists(topParams);

  // Set initial solution and parameters in solution vector
  grpPtr->setX(*startVec);
  *(xVec->getXVec()) = *startVec;
  xVec->getScalar(0) = conParam;

  // Create "C" block in bordered solver
  minusOne = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(1, 1));
  (*minusOne)(0,0) = -1.0;

  // Instantiate bordered solver
  borderedSolver = globalData->locaFactory->createBorderedSolverStrategy(
				   parsedParams,
				   homotopyParams);

  // Determine if underlying group is bordered
  bordered_grp = 
    Teuchos::rcp_dynamic_cast<LOCA::BorderedSystem::AbstractGroup>(grpPtr);
  isBordered = (bordered_grp != Teuchos::null);

  // Create Jacobian operator for bordered solver
  jacOp = Teuchos::rcp(new LOCA::BorderedSolver::JacobianOperator(grpPtr));
}

LOCA::Homotopy::DeflatedGroup::
DeflatedGroup(const LOCA::Homotopy::DeflatedGroup& source,
	      NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    homotopyParams(source.homotopyParams),
    grpPtr(Teuchos::rcp_dynamic_cast<LOCA::Homotopy::AbstractGroup>(source.grpPtr->clone(type))),
    bordered_grp(),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    gradientMultiVec(source.gradientMultiVec, type),
    xVec(),
    fVec(),
    newtonVec(),
    gradientVec(),
    startVec(source.startVec),
    identitySign(source.identitySign),
    solns(source.solns),
    distVec(source.distVec->clone(type)),
    totalDistMultiVec(source.totalDistMultiVec->clone(type)),
    underlyingF(source.underlyingF->clone(type)),
    jacOp(),
    borderedSolver(source.borderedSolver),
    minusOne(Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(*source.minusOne))),
    numSolns(source.numSolns),
    distances(source.distances),
    distProd(source.distProd),
    index_f(1),
    paramVec(source.paramVec),
    conParam(source.conParam),
    conParamID(source.conParamID),
    conParamLabel(source.conParamLabel),
    augmentJacForHomotopyNotImplemented(source.augmentJacForHomotopyNotImplemented),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    isValidGradient(source.isValidGradient),
    isBordered(false)
{
  // Set up multi-vector views
  setupViews();

  // Instantiate bordered solver
  borderedSolver = globalData->locaFactory->createBorderedSolverStrategy(
				   parsedParams,
				   homotopyParams);

  if (type == NOX::ShapeCopy) {
    isValidF = false;
    isValidJacobian = false;
    isValidNewton = false;
    isValidGradient = false;
  }

  // Determine if underlying group is bordered
  bordered_grp = 
    Teuchos::rcp_dynamic_cast<LOCA::BorderedSystem::AbstractGroup>(grpPtr);
  isBordered = (bordered_grp != Teuchos::null);

  // Create Jacobian operator for bordered solver
  jacOp = Teuchos::rcp(new LOCA::BorderedSolver::JacobianOperator(grpPtr));

  // Set blocks in bordered solver
  if (isValidJacobian) {
    borderedSolver->setMatrixBlocksMultiVecConstraint(jacOp, 
						      underlyingF,
						      totalDistMultiVec,
						      minusOne);
    NOX::Abstract::Group::ReturnType status = borderedSolver->initForSolve();
    globalData->locaErrorCheck->checkReturnType(status, 
						"LOCA::Homotopy::DeflatedGroup()");
  }
}


LOCA::Homotopy::DeflatedGroup::
~DeflatedGroup() 
{
}

double
LOCA::Homotopy::DeflatedGroup::
getHomotopyParam() const
{
  return conParam;
}

NOX::Abstract::Group&
LOCA::Homotopy::DeflatedGroup::
operator=(const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Homotopy::DeflatedGroup::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new DeflatedGroup(*this, type));
}

void
LOCA::Homotopy::DeflatedGroup::
setX(const NOX::Abstract::Vector& y)  
{
  const LOCA::MultiContinuation::ExtendedVector& my = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  grpPtr->setX( *(my.getXVec()) );
  *xVec = my;

  resetIsValid();
}

void
LOCA::Homotopy::DeflatedGroup::
computeX(const NOX::Abstract::Group& g, 
	 const NOX::Abstract::Vector& d,
	 double step) 
{
  const LOCA::Homotopy::DeflatedGroup& mg = 
    dynamic_cast<const LOCA::Homotopy::DeflatedGroup&>(g);
  const LOCA::MultiContinuation::ExtendedVector& md = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(d);

  grpPtr->computeX(*(mg.grpPtr), *(md.getXVec()), step);
  xVec->update(1.0, mg.getX(), step, md, 0.0);

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::computeF()";
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

  distProd = 1.0;
  for (int i=0; i<numSolns; i++) {
    distVec->update(1.0, grpPtr->getX(), -1.0, *(solns[i]), 0.0);
    distances[i] = distVec->norm();
    distProd *= distances[i];
  }
  distVec->update(identitySign, grpPtr->getX(), -identitySign, *startVec, 0.0);
  fVec->getXVec()->update(conParam / distProd, grpPtr->getF(), 
			  1.0 - conParam, *distVec, 0.0);
  fVec->getScalar(0) = 0.0;
  (*underlyingF)[0] = grpPtr->getF();

  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  distProd = 1.0;
  totalDistVec->init(0.0);
  for (int i=0; i<numSolns; i++) {
    distVec->update(1.0, grpPtr->getX(), -1.0, *(solns[i]), 0.0);
    distances[i] = distVec->norm();
    distProd *= distances[i];
    totalDistVec->update(-1.0 / (distances[i]*distances[i]), *distVec, 1.0);
  }
  totalDistVec->scale(conParam / distProd);

  NOX::Abstract::Group::ReturnType augHomTest = 
    grpPtr->augmentJacobianForHomotopy(conParam/distProd, identitySign*(1.0-conParam));

  // If it is not implemented, augment the Jacobian during the 
  // applyJacobian() call.
  if (augHomTest == NOX::Abstract::Group::NotDefined)
    augmentJacForHomotopyNotImplemented = true;

  // Set blocks in bordered solver
  borderedSolver->setMatrixBlocksMultiVecConstraint(jacOp, 
						    underlyingF, 
						    totalDistMultiVec,
						    minusOne);
  status = borderedSolver->initForSolve();
  finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::computeGradient()";
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

  // Compute J^T*f for homotopy group
  status = applyJacobianTranspose(*fVec, *gradientVec);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  isValidGradient = true;

  return finalStatus;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
computeNewton(Teuchos::ParameterList& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::computeNewton()";
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
  newtonMultiVec.init(0.0);

  status = applyJacobianInverseMultiVector(params, fMultiVec, 
					   newtonMultiVec);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  newtonMultiVec.scale(-1.0);

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
applyJacobian(const NOX::Abstract::Vector& input,
	      NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobian
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
applyJacobianTranspose(const NOX::Abstract::Vector& input,
		       NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianTranspose
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianTransposeMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
applyJacobianInverse(Teuchos::ParameterList& params, 
		     const NOX::Abstract::Vector& input,
		     NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianInverse
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianInverseMultiVector(params, *mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
applyJacobianMultiVector(const NOX::Abstract::MultiVector& input,
			 NOX::Abstract::MultiVector& result) const 
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::applyJacobianMultiVector()";
  
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
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x = 
    c_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x = 
    c_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    c_result.getScalars();

  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianMultiVector(*input_x, *result_x);

  // If the Jacobian is not augmented for homotopy (i.e. using MFNK)
  // then lets augment it here.
  if (augmentJacForHomotopyNotImplemented)
    result_x->update(1.0-conParam, *input_x, conParam/distProd);

  // Add deflated terms
  if (numSolns > 0) {
    NOX::Abstract::MultiVector::DenseMatrix tmp(1, input.numVectors());
    input_x->multiply(1.0, *totalDistMultiVec, tmp);
    result_x->update(Teuchos::NO_TRANS, 1.0, *underlyingF, tmp, 1.0);
  }

  // Zero out parameter component
  result_param->putScalar(0.0);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
applyJacobianTransposeMultiVector(const NOX::Abstract::MultiVector& input,
				  NOX::Abstract::MultiVector& result) const 
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::applyJacobianTransposeMultiVector()";
  
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
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x = 
    c_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x = 
    c_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    c_result.getScalars();

  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianTransposeMultiVector(*input_x, *result_x);

  // If the Jacobian is not augmented for homotopy (i.e. using MFNK)
  // then lets augment it here.
  if (augmentJacForHomotopyNotImplemented)
    result_x->update(1.0-conParam, *input_x, conParam/distProd);

  // Add deflated terms
  if (numSolns > 0) {
    NOX::Abstract::MultiVector::DenseMatrix tmp(1, input.numVectors());
    input_x->multiply(1.0, *underlyingF, tmp);
    result_x->update(Teuchos::NO_TRANS, 1.0, *totalDistMultiVec, tmp, 1.0);
  }

  // Zero out parameter component
  result_param->putScalar(0.0);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
applyJacobianInverseMultiVector(Teuchos::ParameterList& params,
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const 
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::applyJacobianInverseMultiVector()";
  
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
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x = 
    c_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x = 
    c_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    c_result.getScalars();

  NOX::Abstract::Group::ReturnType status;
  if (numSolns > 0) {
    // Call bordered solver applyInverse method
    status = 
      borderedSolver->applyInverse(params, input_x.get(), input_param.get(), 
				   *result_x, *result_param);
  }
  else {
    status = 
      grpPtr->applyJacobianInverseMultiVector(params, *input_x, *result_x);
    result_param->putScalar(0.0);
  }

  return status;
}

bool
LOCA::Homotopy::DeflatedGroup::
isF() const 
{
  return isValidF;
}

bool
LOCA::Homotopy::DeflatedGroup::
isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Homotopy::DeflatedGroup::
isGradient() const 
{
  return isValidGradient;
}

bool
LOCA::Homotopy::DeflatedGroup::
isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Homotopy::DeflatedGroup::
getX() const 
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::Homotopy::DeflatedGroup::
getF() const 
{
  return *fVec;
}

double
LOCA::Homotopy::DeflatedGroup::
getNormF() const 
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::Homotopy::DeflatedGroup::
getGradient() const 
{
  return *gradientVec;
}

const NOX::Abstract::Vector&
LOCA::Homotopy::DeflatedGroup::
getNewton() const 
{
  return *newtonVec;
}
  
Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::
getXPtr() const 
{
  return xVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::
getFPtr() const 
{
  return fVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::
getGradientPtr() const 
{
  return gradientVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::
getNewtonPtr() const 
{
  return newtonVec;
}

double
LOCA::Homotopy::DeflatedGroup::
getNormNewtonSolveResidual() const 
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::MultiContinuation::ExtendedVector residual = *fVec;
  
  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual = residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Homotopy::DeflatedGroup::
getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Homotopy::DeflatedGroup::
getUnderlyingGroup()
{
  return grpPtr;
}

void
LOCA::Homotopy::DeflatedGroup::
copy(const NOX::Abstract::Group& src) 
{

  const LOCA::Homotopy::DeflatedGroup& source = 
    dynamic_cast<const LOCA::Homotopy::DeflatedGroup&>(src);

  // Protect against A = A
  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    homotopyParams = source.homotopyParams;
    grpPtr->copy(*source.grpPtr);
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    gradientMultiVec = source.gradientMultiVec;
    startVec = source.startVec;
    identitySign = source.identitySign;
    solns = source.solns;
    *distVec = *source.distVec;
    *totalDistMultiVec = *source.totalDistMultiVec;
    *underlyingF = *source.underlyingF;
    numSolns = source.numSolns;
    distances = source.distances;
    distProd = source.distProd;
    index_f = source.index_f;
    paramVec = source.paramVec;
    conParam = source.conParam;
    conParamID = source.conParamID;
    augmentJacForHomotopyNotImplemented = source.augmentJacForHomotopyNotImplemented;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;

    // set up views again just to be safe
    setupViews();

    // Instantiate bordered solver
    borderedSolver = 
      globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
							    homotopyParams);

    // Set blocks in bordered solver
    if (isValidJacobian) {
      borderedSolver->setMatrixBlocksMultiVecConstraint(jacOp, 
							underlyingF,
							totalDistMultiVec,
							minusOne);
      NOX::Abstract::Group::ReturnType status = borderedSolver->initForSolve();
      globalData->locaErrorCheck->checkReturnType(status, 
						  "LOCA::Homotopy::copy()");
    }
  }
}

void
LOCA::Homotopy::DeflatedGroup::
setParamsMulti(const std::vector<int>& paramIDs, 
	       const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);

  for (unsigned int i=0; i<paramIDs.size(); i++) {
    paramVec[paramIDs[i]] = vals(i,0);
    if (paramIDs[i] == conParamID)
      conParam = vals(i,0);
  }

  resetIsValid();
}

void
LOCA::Homotopy::DeflatedGroup::
setParams(const LOCA::ParameterVector& p)
{
  grpPtr->setParams(p);
  xVec->getScalar(0) = p[conParamID];

  resetIsValid();
}

void
LOCA::Homotopy::DeflatedGroup::
setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
  paramVec[paramID] = val;
  if (paramID == conParamID)
    conParam = val;

  resetIsValid();
}

void
LOCA::Homotopy::DeflatedGroup::
setParam(std::string paramID, double val)
{
  int id = paramVec.getIndex(paramID);
  setParam(id, val);
}

const LOCA::ParameterVector&
LOCA::Homotopy::DeflatedGroup::
getParams() const
{
  return paramVec;
}

double
LOCA::Homotopy::DeflatedGroup::
getParam(int paramID) const{
  return paramVec[paramID];
}

double
LOCA::Homotopy::DeflatedGroup::
getParam(std::string paramID) const
{
  return paramVec.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::
computeDfDpMulti(const std::vector<int>& paramIDs, 
		 NOX::Abstract::MultiVector& dfdp, 
		 bool isValid_F)
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::computeDfDpMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast dfdp to an extended multi-vec
  LOCA::MultiContinuation::ExtendedMultiVector& e_dfdp = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(dfdp);

  Teuchos::RCP<NOX::Abstract::MultiVector> dfdp_x = 
    e_dfdp.getXMultiVec();

  if (!isValid_F) {
    status = grpPtr->computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    (*dfdp_x)[0].update(conParam / distProd, grpPtr->getF(), 0.0);
  }

  std::vector<int> index_c, index_p, p_IDs;
  index_p.push_back(0);
  for (unsigned int i=0; i<paramIDs.size(); i++) {
    if (paramIDs[i] == conParamID)
      index_c.push_back(i+1);
    else {
      index_p.push_back(i+1);
      p_IDs.push_back(paramIDs[i]);
    }
  }

  Teuchos::RCP<NOX::Abstract::MultiVector> dfdp_c, dfdp_p;
  if (index_c.size() > 0) {
    dfdp_c = dfdp_x->subView(index_c);
    distVec->update(1.0, grpPtr->getX(), -1.0, *startVec, 0.0);
    
    (*dfdp_c)[0].update(1.0/distProd, grpPtr->getF(), -1.0, *distVec, 0.0);
  }

  if (index_p.size() > 1) {
    dfdp_p = dfdp_x->subView(index_p);
    status = grpPtr->computeDfDpMulti(p_IDs, *(dfdp_p), true);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    for (unsigned int i=0; i<p_IDs.size(); i++) {
      (*dfdp_p)[i+1].scale(conParam / distProd);
    }
  }

  // Set param component to 0
  e_dfdp.getScalars()->putScalar(0.0);

  return finalStatus;
}

void
LOCA::Homotopy::DeflatedGroup::
preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->preProcessContinuationStep(stepStatus);
}

void
LOCA::Homotopy::DeflatedGroup::
postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->postProcessContinuationStep(stepStatus);
}

void
LOCA::Homotopy::DeflatedGroup::
projectToDraw(const NOX::Abstract::Vector& x,
	      double *px) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);

  grpPtr->projectToDraw(*mx.getXVec(), px);
  px[grpPtr->projectToDrawDimension()] = mx.getScalar(0);
}

int
LOCA::Homotopy::DeflatedGroup::
projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension() + 1;
}

double
LOCA::Homotopy::DeflatedGroup::
computeScaledDotProduct(const NOX::Abstract::Vector& a,
			const NOX::Abstract::Vector& b) const
{
  const LOCA::MultiContinuation::ExtendedVector& ma = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(a);
  const LOCA::MultiContinuation::ExtendedVector& mb = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(b);

  double val = grpPtr->computeScaledDotProduct(*ma.getXVec(), *mb.getXVec());
  val += ma.getScalar(0) * mb.getScalar(0);

  return val;
}

void
LOCA::Homotopy::DeflatedGroup::
printSolution(const double conParam_) const
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for homotopy parameter = " << 
      globalData->locaUtils->sciformat(conParam_) << std::endl;
  }
  grpPtr->printSolution(conParam_);
  return;
}

void
LOCA::Homotopy::DeflatedGroup::
printSolution(const NOX::Abstract::Vector& x,
	      const double conParam_) const
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for homotopy parameter = " << 
      globalData->locaUtils->sciformat(conParam_) << std::endl;
  }
  grpPtr->printSolution(x, conParam_);
  return;
}

void
LOCA::Homotopy::DeflatedGroup::
scaleVector(NOX::Abstract::Vector& x) const
{
  LOCA::MultiContinuation::ExtendedVector& mx = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector&>(x);

  grpPtr->scaleVector(*mx.getXVec());
}

int
LOCA::Homotopy::DeflatedGroup::
getBorderedWidth() const
{
  int my_width = 1;
  if (isBordered)
    return my_width + bordered_grp->getBorderedWidth();
  else
    return my_width;
}

Teuchos::RCP<const NOX::Abstract::Group>
LOCA::Homotopy::DeflatedGroup::
getUnborderedGroup() const
{
  if (isBordered)
    return bordered_grp->getUnborderedGroup();
  else
    return grpPtr;
}

bool
LOCA::Homotopy::DeflatedGroup::
isCombinedAZero() const
{
  return false;  // A is always considered non-zero 
}

bool
LOCA::Homotopy::DeflatedGroup::
isCombinedBZero() const
{
  return false; // B is always considered non-zero 
}

bool
LOCA::Homotopy::DeflatedGroup::
isCombinedCZero() const
{
  return false;  // C is always considered non-zero 
}

void
LOCA::Homotopy::DeflatedGroup::
extractSolutionComponent(const NOX::Abstract::MultiVector& v,
			 NOX::Abstract::MultiVector& v_x) const
{
  // cast v to an extended multi-vec
  const LOCA::MultiContinuation::ExtendedMultiVector& mc_v = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(v);
  
  // get solution component
  Teuchos::RCP<const NOX::Abstract::MultiVector> mc_v_x =
    mc_v.getXMultiVec();

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    v_x = *mc_v_x;
    return;
  }

  // Extract solution component from mc_v_x
  bordered_grp->extractSolutionComponent(*mc_v_x, v_x);
}

void
LOCA::Homotopy::DeflatedGroup::
extractParameterComponent(bool use_transpose,
			  const NOX::Abstract::MultiVector& v,
			  NOX::Abstract::MultiVector::DenseMatrix& v_p) const
{
  // cast v to an extended multi-vec
  const LOCA::MultiContinuation::ExtendedMultiVector& mc_v = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(v);
  
  // get solution and parameter components
  Teuchos::RCP<const NOX::Abstract::MultiVector> mc_v_x =
    mc_v.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> mc_v_p =
    mc_v.getScalars();

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    if (!use_transpose)
      v_p.assign(*mc_v_p);
    else
      for (int j=0; j<v_p.numCols(); j++)
	for (int i=0; i<v_p.numRows(); i++)
	  v_p(i,j) = (*mc_v_p)(j,i);
    return;
  }

  int w = bordered_grp->getBorderedWidth();
  if (!use_transpose) {
    // Split v_p into 2 block rows, the top to store mc_v_x_p and the bottom
    // to store mc_v_p
    int num_cols = v_p.numCols();
    NOX::Abstract::MultiVector::DenseMatrix v_p_1(Teuchos::View, v_p,
						  w, num_cols, 0, 0);
    NOX::Abstract::MultiVector::DenseMatrix v_p_2(Teuchos::View, v_p,
						  1, num_cols, w, 0);

    // Decompose mc_v_x
    bordered_grp->extractParameterComponent(use_transpose,*mc_v_x, v_p_1);
    v_p_2.assign(*mc_v_p);
  }
  else {
    // Split v_p into 2 block columns, the first to store mc_v_x_p^t and the 
    // the second to store mc_v_p^T
    int num_rows = v_p.numRows();
    NOX::Abstract::MultiVector::DenseMatrix v_p_1(Teuchos::View, v_p,
						  num_rows, w, 0, 0);
    NOX::Abstract::MultiVector::DenseMatrix v_p_2(Teuchos::View, v_p,
						  num_rows, 1, 0, w);

    // Decompose mc_v_x
    bordered_grp->extractParameterComponent(use_transpose,*mc_v_x, v_p_1);
    for (int j=0; j<1; j++)
      for (int i=0; i<num_rows; i++)
	v_p_2(i,j) = (*mc_v_p)(j,i);
  }
}

void
LOCA::Homotopy::DeflatedGroup::
loadNestedComponents(const NOX::Abstract::MultiVector& v_x,
		     const NOX::Abstract::MultiVector::DenseMatrix& v_p,
		     NOX::Abstract::MultiVector& v) const
{
  // cast X to an extended multi-vec
  LOCA::MultiContinuation::ExtendedMultiVector& mc_v = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector&>(v);

  // get solution and parameter components
  Teuchos::RCP<NOX::Abstract::MultiVector> mc_v_x =
    mc_v.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> mc_v_p =
    mc_v.getScalars();

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    *mc_v_x = v_x;
    mc_v_p->assign(v_p);
    return;
  }

  // split v_p
  int num_cols = v_p.numCols();
  int w = bordered_grp->getBorderedWidth();
  NOX::Abstract::MultiVector::DenseMatrix v_p_1(Teuchos::View, v_p,
						w, num_cols, 0, 0);
  NOX::Abstract::MultiVector::DenseMatrix v_p_2(Teuchos::View, v_p,
						1, num_cols, w, 0);

  // load v_x, v_p_1 into mc_v_x
  bordered_grp->loadNestedComponents(v_x, v_p_1, *mc_v_x);

  // load v_p_2 into mc_v_p
  mc_v_p->assign(v_p_2);
}

void
LOCA::Homotopy::DeflatedGroup::
fillA(NOX::Abstract::MultiVector& A) const
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::fillA";

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_A = 
    underlyingF;

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    A = *my_A;
    return;
  }

  // Create views for underlying group
  int w = bordered_grp->getBorderedWidth();
  std::vector<int> idx1(w);
  for (int i=0; i<w; i++)
    idx1[i] = i;
  Teuchos::RCP<NOX::Abstract::MultiVector> underlyingA = 
    A.subView(idx1);

  // Fill A block in underlying group
  bordered_grp->fillA(*underlyingA);

  // Create views for my blocks
  std::vector<int> idx2(1);
  for (int i=0; i<1; i++)
    idx2[i] = w+i;
  Teuchos::RCP<NOX::Abstract::MultiVector> my_A_x = 
    A.subView(idx2);

  // Extract solution component from my_A and store in A
  bordered_grp->extractSolutionComponent(*my_A, *my_A_x);
}

void
LOCA::Homotopy::DeflatedGroup::
fillB(NOX::Abstract::MultiVector& B) const
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::fillB";

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_B =
    totalDistMultiVec;

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    B = *my_B;
    return;
  }

  // Create views for underlying group
  int w = bordered_grp->getBorderedWidth();
  std::vector<int> idx1(w);
  for (int i=0; i<w; i++)
    idx1[i] = i;
  Teuchos::RCP<NOX::Abstract::MultiVector> underlyingB = 
    B.subView(idx1);

  // Combine blocks in underlying group
  bordered_grp->fillB(*underlyingB);

  // Create views for my blocks
  std::vector<int> idx2(2);
  for (int i=0; i<1; i++)
    idx2[i] = w+i;
  Teuchos::RCP<NOX::Abstract::MultiVector> my_B_x = 
    B.subView(idx2);

  // Extract solution component from my_B and store in B
  bordered_grp->extractSolutionComponent(*my_B, *my_B_x);
}

void
LOCA::Homotopy::DeflatedGroup::
fillC(NOX::Abstract::MultiVector::DenseMatrix& C) const
{
  std::string callingFunction = 
    "LOCA::Homotopy::DeflatedGroup::fillC";

  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> my_C = 
    minusOne;

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    C.assign(*my_C);
    return;
  }

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_B = 
    totalDistMultiVec;

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_A = 
    underlyingF;
  
  // Create views for underlying group
  int w = bordered_grp->getBorderedWidth();
  NOX::Abstract::MultiVector::DenseMatrix underlyingC(Teuchos::View, C,
						      w, w, 0, 0);

  // Combine blocks in underlying group
  bordered_grp->fillC(underlyingC);

  // Create views for my blocks
  NOX::Abstract::MultiVector::DenseMatrix my_A_p(Teuchos::View, C,
						 w, 1, 0, w);
  NOX::Abstract::MultiVector::DenseMatrix my_B_p(Teuchos::View, C,
						 1, w, w, 0);
  NOX::Abstract::MultiVector::DenseMatrix my_CC(Teuchos::View, C,
						1, 1, w, w);

  // Extract solution component from my_A and store in my_A_p
  bordered_grp->extractParameterComponent(false, *my_A, my_A_p);

  // Extract solution component from my_B and store in my_B_p
  bordered_grp->extractParameterComponent(true, *my_B, my_B_p);

  // Copy in my_C
  my_CC.assign(*my_C);
}

void
LOCA::Homotopy::DeflatedGroup::
resetIsValid() 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
}

void
LOCA::Homotopy::DeflatedGroup::
setupViews()
{
  index_f[0] = 0;

  xVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(xMultiVec.getVector(0),true);
  fVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(fMultiVec.getVector(0),true);
  newtonVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(newtonMultiVec.getVector(0),true);
  gradientVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(gradientMultiVec.getVector(0),true);
  totalDistVec = Teuchos::rcp(&(*totalDistMultiVec)[0], false);
}

void
LOCA::Homotopy::DeflatedGroup::
setHomotopyParam(double val) 
{
  xVec->getScalar(0) = val;
  paramVec[conParamID] = val;

  resetIsValid();
}

void
LOCA::Homotopy::DeflatedGroup::
setStepperParameters(Teuchos::ParameterList& topParams)
{
  Teuchos::ParameterList& params = topParams.sublist("LOCA");
  
  // Create the stepper sublist
  Teuchos::ParameterList& stepperList = params.sublist("Stepper");
  stepperList.set("Continuation Method", "Natural");
  stepperList.set("Continuation Parameter", conParamLabel);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  stepperList.set("Min Value", -1.0);
  stepperList.set("Max Steps", 50);

  // Create predictor sublist
  Teuchos::ParameterList& predictorList = params.sublist("Predictor");
  predictorList.set("Method", "Constant");

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = params.sublist("Step Size");
  //stepSizeList.set("Method", "Constant");
  stepSizeList.set("Method", "Adaptive");
  stepSizeList.set("Initial Step Size", 0.1);
  stepSizeList.set("Min Step Size", 1.0e-2);
  stepSizeList.set("Max Step Size", 1.0);
  stepSizeList.set("Aggressiveness", 0.5);
  return;
}
