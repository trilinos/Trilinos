// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_BorderedSolver_JacobianOperator.H"

LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup(
       const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& conParams,
       const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& g,
       const Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>& constraints,
       const std::vector<int>& paramIDs,
       bool skip_dfdp)
  : LOCA::Abstract::Group(global_data),
    globalData(global_data),
    parsedParams(topParams),
    constraintParams(conParams),
    grpPtr(g),
    bordered_grp(),
    constraintsPtr(constraints),
    numParams(paramIDs.size()),
    xMultiVec(globalData, g->getX(), 1, numParams, NOX::DeepCopy),
    fMultiVec(globalData, g->getX(), numParams+1, numParams, NOX::ShapeCopy),
    newtonMultiVec(globalData, g->getX(), 1, numParams, NOX::ShapeCopy),
    gradientMultiVec(globalData, g->getX(), 1, numParams, NOX::ShapeCopy),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    gradientVec(),
    jacOp(),
    borderedSolver(),
    index_f(1),
    index_dfdp(numParams),
    constraintParamIDs(paramIDs),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false),
    isValidGradient(false),
    isBordered(false),
    skipDfDp(skip_dfdp)
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
  borderedSolver = globalData->locaFactory->createBorderedSolverStrategy(
                   parsedParams,
                   constraintParams);

  // Determine if underlying group is bordered
  bordered_grp =
    Teuchos::rcp_dynamic_cast<LOCA::BorderedSystem::AbstractGroup>(grpPtr);
  isBordered = (bordered_grp != Teuchos::null);

  // Create Jacobian operator for bordered solver
  jacOp = Teuchos::rcp(new LOCA::BorderedSolver::JacobianOperator(grpPtr));
}

LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup(
             const LOCA::MultiContinuation::ConstrainedGroup& source,
             NOX::CopyType type)
  : LOCA::Abstract::Group(source),
    globalData(source.globalData),
    parsedParams(source.parsedParams),
    constraintParams(source.constraintParams),
    grpPtr(Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::AbstractGroup>(source.grpPtr->clone(type))),
    bordered_grp(),
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
    jacOp(),
    borderedSolver(source.borderedSolver),
    index_f(1),
    index_dfdp(numParams),
    constraintParamIDs(source.constraintParamIDs),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    isValidGradient(source.isValidGradient),
    isBordered(false),
    skipDfDp(source.skipDfDp)
{
  // Set up multi-vector views
  setupViews();

  // Instantiate bordered solver
  borderedSolver = globalData->locaFactory->createBorderedSolverStrategy(
                   parsedParams,
                   constraintParams);

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
    if (skipDfDp)
      borderedSolver->setMatrixBlocks(jacOp,
                      Teuchos::null,
                      constraintsPtr,
                      dfdpMultiVec->getScalars());
    else
      borderedSolver->setMatrixBlocks(jacOp,
                      dfdpMultiVec->getXMultiVec(),
                      constraintsPtr,
                      dfdpMultiVec->getScalars());
    NOX::Abstract::Group::ReturnType status = borderedSolver->initForSolve();
    globalData->locaErrorCheck->checkReturnType(status,
                        "LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup()");
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
  constraintsPtr->setParam(constraintParamIDs[i],val);

  resetIsValid();
}

double
LOCA::MultiContinuation::ConstrainedGroup::getConstraintParameter(int i) const
{
  return grpPtr->getParam(constraintParamIDs[i]);
}

const std::vector<int>&
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

Teuchos::RCP<NOX::Abstract::Group>
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

  std::string callingFunction =
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
  if (isValidJacobian && grpPtr->isJacobian())
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::MultiContinuation::ConstrainedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  if (!skipDfDp) {
    status = grpPtr->computeDfDpMulti(constraintParamIDs,
                      *fMultiVec.getXMultiVec(),
                      isValidF);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // We need to compute the constraint derivatives before computing the
  // Jacobian, because the constraint derivatives might involve derivatives
  // of the Jacobian, and finite differencing may invalidate the Jacobian

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

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Set blocks in bordered solver
  if (skipDfDp)
    borderedSolver->setMatrixBlocks(jacOp,
                    Teuchos::null,
                    constraintsPtr,
                    dfdpMultiVec->getScalars());
  else
    borderedSolver->setMatrixBlocks(jacOp,
                    dfdpMultiVec->getXMultiVec(),
                    constraintsPtr,
                    dfdpMultiVec->getScalars());
  status = borderedSolver->initForSolve();
  finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeGradient()
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
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
                           Teuchos::ParameterList& params)
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
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
  newtonMultiVec.init(0.0);

  status = applyJacobianInverseMultiVector(params, *ffMultiVec,
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
LOCA::MultiContinuation::ConstrainedGroup::applyJacobian(
                      const NOX::Abstract::Vector& input,
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
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTranspose(
                      const NOX::Abstract::Vector& input,
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
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverse(
                      Teuchos::ParameterList& params,
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
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianMultiVector(
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
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
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x =
    c_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x =
    c_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param =
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
  std::string callingFunction =
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
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x =
    c_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = c_input.getScalars();

  // Get references to x, param components of result vector
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x =
    c_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param =
    c_result.getScalars();

  // Call bordered solver applyTranspose method
  NOX::Abstract::Group::ReturnType status =
    borderedSolver->applyTranspose(*input_x, *input_param, *result_x,
                   *result_param);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseMultiVector(
                     Teuchos::ParameterList& params,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

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

  // Call bordered solver applyInverse method
  status = borderedSolver->initForSolve();
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);
  status =
    borderedSolver->applyInverse(params, input_x.get(), input_param.get(),
                 *result_x, *result_param);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

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

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getXPtr() const
{
  return xVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getFPtr() const
{
  return fVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getGradientPtr() const
{
  return gradientVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getNewtonPtr() const
{
  return newtonVec;
}

double
LOCA::MultiContinuation::ConstrainedGroup::getNormNewtonSolveResidual() const
{
  std::string callingFunction =
    "LOCA::MultiContinuation::ConstrainedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::MultiContinuation::ExtendedVector residual = *fVec;

  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual = residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup()
{
  return grpPtr;
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
    skipDfDp = source.skipDfDp;

    // set up views again just to be safe
    setupViews();

    // Instantiate bordered solver
    borderedSolver =
      globalData->locaFactory->createBorderedSolverStrategy(
                   parsedParams,
                   constraintParams);

    // Set blocks in bordered solver
    if (isValidJacobian) {
      if (skipDfDp)
    borderedSolver->setMatrixBlocks(jacOp,
                    Teuchos::null,
                    constraintsPtr,
                    dfdpMultiVec->getScalars());
      else
    borderedSolver->setMatrixBlocks(jacOp,
                    dfdpMultiVec->getXMultiVec(),
                    constraintsPtr,
                    dfdpMultiVec->getScalars());
      NOX::Abstract::Group::ReturnType status = borderedSolver->initForSolve();
      globalData->locaErrorCheck->checkReturnType(status,
                          "LOCA::MultiContinuation::ConstrainedGroup::copy()");
    }
  }
}

void
LOCA::MultiContinuation::ConstrainedGroup::setParamsMulti(
             const std::vector<int>& paramIDs,
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
LOCA::MultiContinuation::ConstrainedGroup::setParam(std::string paramID, double val)
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
LOCA::MultiContinuation::ConstrainedGroup::getParam(std::string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeDfDpMulti(
                         const std::vector<int>& paramIDs,
                         NOX::Abstract::MultiVector& dfdp,
                         bool isValid_F)
{
  std::string callingFunction =
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
LOCA::MultiContinuation::ConstrainedGroup::preProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->preProcessContinuationStep(stepStatus);
  constraintsPtr->preProcessContinuationStep(stepStatus);
}

void
LOCA::MultiContinuation::ConstrainedGroup::postProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->postProcessContinuationStep(stepStatus);
  constraintsPtr->postProcessContinuationStep(stepStatus);
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
    const LOCA::ParameterVector& pvec = grpPtr->getParams();
    globalData->locaUtils->out() << "\tPrinting constraint parameters:\n";
    for (int i=0; i<numParams; i++)
      globalData->locaUtils->out() << "\t\t" <<
    pvec.getLabel(constraintParamIDs[i]) << " = " <<
    globalData->locaUtils->sciformat(mx.getScalar(i)) << std::endl;
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

int
LOCA::MultiContinuation::ConstrainedGroup::getBorderedWidth() const
{
  int my_width = numParams;
  if (isBordered)
    return my_width + bordered_grp->getBorderedWidth();
  else
    return my_width;
}

Teuchos::RCP<const NOX::Abstract::Group>
LOCA::MultiContinuation::ConstrainedGroup::getUnborderedGroup() const
{
  if (isBordered)
    return bordered_grp->getUnborderedGroup();
  else
    return grpPtr;
}

bool
LOCA::MultiContinuation::ConstrainedGroup::isCombinedAZero() const
{
  return false;  // A is always considered non-zero (A = df/dp)
}

bool
LOCA::MultiContinuation::ConstrainedGroup::isCombinedBZero() const
{
  if (isBordered)
    return constraintsPtr->isDXZero() && bordered_grp->isCombinedBZero();
  else
    return constraintsPtr->isDXZero();
}

bool
LOCA::MultiContinuation::ConstrainedGroup::isCombinedCZero() const
{
  return false;  // C is always considered non-zero (C = dg/dp)
}

// void
// LOCA::MultiContinuation::ConstrainedGroup::decomposeMultiVec(
//                bool use_transpose,
//                            const NOX::Abstract::MultiVector& v,
//                            NOX::Abstract::MultiVector& v_x,
//                            NOX::Abstract::MultiVector::DenseMatrix& v_p) const
// {
//   // cast v to an extended multi-vec
//   const LOCA::MultiContinuation::ExtendedMultiVec& mc_v =
//     dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVec&>(v);

//   // break mc_v into solution and parameter components
//   Teuchos::RCP<const NOX::Abstract::MultiVector> mc_v_x =
//     mc_v.getXMultiVec();
//   Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> mc_v_p =
//     mc_v.getScalars();

//   // If the underlying system isn't bordered, we're done
//   if (!isBordered) {
//     v_x = *mc_v_x;
//     if (!use_transpose)
//       v_p.assign(*mc_v_p);
//     else
//       for (int j=0; j<v_p.numCols(); j++)
//     for (int i=0; i<v_p.numRows(); i++)
//       v_p(i,j) = (*mc_v_p)(j,i);
//     return;
//   }

//   int w = bordered_grp->getBorderedWidth();
//   if (!use_transpose) {
//     // Split v_p into 2 block rows, the top to store mc_v_x_p and the bottom
//     // to store mc_v_p

//     int num_cols = v_p.numCols();
//     NOX::Abstract::MultiVector::DenseMatrix v_p_1(Teuchos::View, v_p,
//                           w, num_cols, 0, 0);
//     NOX::Abstract::MultiVector::DenseMatrix v_p_2(Teuchos::View, v_p,
//                           numParams, num_cols, w, 0);

//     // Decompose mc_v_x
//     bordered_grp->decomposeMultiVec(use_transpose,*mc_v_x, v_x, v_p_1);
//     v_p_2.assign(*mc_v_p);
//   }
//   else {
//     // Split v_p into 2 block columns, the first to store mc_v_x_p^t and the
//     // the second to store mc_v_p^T
//     int num_rows = v_p.numRows();
//     NOX::Abstract::MultiVector::DenseMatrix v_p_1(Teuchos::View, v_p,
//                           num_rows, w, 0, 0);
//     NOX::Abstract::MultiVector::DenseMatrix v_p_2(Teuchos::View, v_p,
//                           num_rows, numParams, 0, w);

//     // Decompose mc_v_x
//     bordered_grp->decomposeMultiVec(use_transpose,*mc_v_x, v_x, v_p_1);
//     for (int j=0; j<numParams; j++)
//       for (int i=0; i<num_rows; i++)
//     v_p_2(i,j) = (*mc_v_p)(j,i);
//   }
// }

// void
// LOCA::MultiContinuation::ConstrainedGroup::combineBlocks(
//                          NOX::Abstract::MultiVector& A,
//                  NOX::Abstract::MultiVector& B,
//                          NOX::Abstract::MultiVector::DenseMatrix& C) const
// {
//   std::string callingFunction =
//     "LOCA::MultiContinuation::ConstrainedGroup::combineBlocks";

//   Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> constraints_mvdx = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(constraints);
//   if (constraints_mvdx == Teuchos::null)
//     global_data.locaErrorCheck->throwError(
//                 callingFunction,
//                 std::string("Constraints object must be of type") +
//                 std::string("ConstraintInterfaceMVDX"));

//   Teuchos::RCP<const NOX::Abstract::MultiVector> my_A =
//     dfdpMultiVec->getXMultiVec();
//   Teuchos::RCP<const NOX::Abstract::MultiVector> my_B =
//     Teuchos::rcp(constraints_mvdx->getDX(),false);
//   Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> my_C =
//     dfdpMultiVec->getScalars();

//   // If the underlying system isn't bordered, we're done
//   if (!isBordered) {
//     A = *my_A;
//     B = *my_B;
//     C = *my_C;
//     return;
//   }

//   // Create views for underlying group
//   int w = bordered_grp->getBorderedWidth();
//   std::vector<int> idx1(w);
//   for (int i=0; i<w; i++)
//     idx1[i] = i;
//   Teuchos::RCP<NOX::Abstract::MultiVector> underlyingA =
//     A.subView(idx1);
//   Teuchos::RCP<NOX::Abstract::MultiVector> underlyingB =
//     B.subView(idx1);
//   NOX::Abstract::MultiVector::DenseMatrix underlyingC(Teuchos::View, C,
//                               w, w, 0, 0);

//   // Combine blocks in underlying group
//   bordered_grp->combineBlocks(underlyingA, underlyingB, underlyingC);

//   // Create views for my blocks
//   std::vector<int> idx2(numParams);
//   for (int i=0; i<numParams; i++)
//     idx2[i] = w+i;
//   Teuchos::RCP<NOX::Abstract::MultiVector> my_A_x =
//     A.subView(idx2);
//   Teuchos::RCP<NOX::Abstract::MultiVector> my_B_x =
//     B.subView(idx2);
//   NOX::Abstract::MultiVector::DenseMatrix my_A_p(Teuchos::View, C,
//                          w, numParams, 0, w);
//   NOX::Abstract::MultiVector::DenseMatrix my_B_p(Teuchos::View, C,
//                          numParams, w, w, 0);
//   NOX::Abstract::MultiVector::DenseMatrix my_CC(Teuchos::View, C,
//                         numParams, numParams, w, w);

//   // Split my A into solution and parameter components
//   bordered_grp->decomposeMultiVec(false, my_A, my_A_x, my_A_p);

//   // Split my B into solution and parameter components
//   bordered_grp->decomposeMultiVec(true, my_B, my_B_x, my_B_p);

//   // Copy in my_C
//   my_CC.assign(*my_C);
// }

void
LOCA::MultiContinuation::ConstrainedGroup::extractSolutionComponent(
                            const NOX::Abstract::MultiVector& v,
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
LOCA::MultiContinuation::ConstrainedGroup::extractParameterComponent(
               bool use_transpose,
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
                          numParams, num_cols, w, 0);

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
                          num_rows, numParams, 0, w);

    // Decompose mc_v_x
    bordered_grp->extractParameterComponent(use_transpose,*mc_v_x, v_p_1);
    for (int j=0; j<numParams; j++)
      for (int i=0; i<num_rows; i++)
    v_p_2(i,j) = (*mc_v_p)(j,i);
  }
}

void
LOCA::MultiContinuation::ConstrainedGroup::loadNestedComponents(
               const NOX::Abstract::MultiVector& v_x,
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
                        numParams, num_cols, w, 0);

  // load v_x, v_p_1 into mc_v_x
  bordered_grp->loadNestedComponents(v_x, v_p_1, *mc_v_x);

  // load v_p_2 into mc_v_p
  mc_v_p->assign(v_p_2);
}

void
LOCA::MultiContinuation::ConstrainedGroup::fillA(
                                     NOX::Abstract::MultiVector& A) const
{
  std::string callingFunction =
    "LOCA::MultiContinuation::ConstrainedGroup::fillA";

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_A =
    dfdpMultiVec->getXMultiVec();

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
  std::vector<int> idx2(numParams);
  for (int i=0; i<numParams; i++)
    idx2[i] = w+i;
  Teuchos::RCP<NOX::Abstract::MultiVector> my_A_x =
    A.subView(idx2);

  // Extract solution component from my_A and store in A
  bordered_grp->extractSolutionComponent(*my_A, *my_A_x);
}

void
LOCA::MultiContinuation::ConstrainedGroup::fillB(
                                      NOX::Abstract::MultiVector& B) const
{
  std::string callingFunction =
    "LOCA::MultiContinuation::ConstrainedGroup::fillB";

  bool isZeroB = constraintsPtr->isDXZero();
  Teuchos::RCP<const NOX::Abstract::MultiVector> my_B;

  if (!isZeroB) {
    Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> constraints_mvdx = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(constraintsPtr);
    if (constraints_mvdx == Teuchos::null)
      globalData->locaErrorCheck->throwError(
                callingFunction,
                std::string("Constraints object must be of type") +
                std::string("ConstraintInterfaceMVDX"));

    my_B = Teuchos::rcp(constraints_mvdx->getDX(),false);
  }

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    if (isZeroB)
      B.init(0.0);
    else
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
  std::vector<int> idx2(numParams);
  for (int i=0; i<numParams; i++)
    idx2[i] = w+i;
  Teuchos::RCP<NOX::Abstract::MultiVector> my_B_x =
    B.subView(idx2);

  // Extract solution component from my_B and store in B
  if (isZeroB)
    my_B_x->init(0.0);
  else
    bordered_grp->extractSolutionComponent(*my_B, *my_B_x);
}

void
LOCA::MultiContinuation::ConstrainedGroup::fillC(
                         NOX::Abstract::MultiVector::DenseMatrix& C) const
{
  std::string callingFunction =
    "LOCA::MultiContinuation::ConstrainedGroup::fillC";

  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> my_C =
    dfdpMultiVec->getScalars();

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    C.assign(*my_C);
    return;
  }

  bool isZeroB = constraintsPtr->isDXZero();
  Teuchos::RCP<const NOX::Abstract::MultiVector> my_B;

  if (!isZeroB) {
    Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> constraints_mvdx = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(constraintsPtr);
    if (constraints_mvdx == Teuchos::null)
      globalData->locaErrorCheck->throwError(
                callingFunction,
                std::string("Constraints object must be of type") +
                std::string("ConstraintInterfaceMVDX"));

    my_B = Teuchos::rcp(constraints_mvdx->getDX(),false);
  }

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_A =
    dfdpMultiVec->getXMultiVec();

  // Create views for underlying group
  int w = bordered_grp->getBorderedWidth();
  NOX::Abstract::MultiVector::DenseMatrix underlyingC(Teuchos::View, C,
                              w, w, 0, 0);

  // Combine blocks in underlying group
  bordered_grp->fillC(underlyingC);

  // Create views for my blocks
  NOX::Abstract::MultiVector::DenseMatrix my_A_p(Teuchos::View, C,
                         w, numParams, 0, w);
  NOX::Abstract::MultiVector::DenseMatrix my_B_p(Teuchos::View, C,
                         numParams, w, w, 0);
  NOX::Abstract::MultiVector::DenseMatrix my_CC(Teuchos::View, C,
                        numParams, numParams, w, w);

  // Extract solution component from my_A and store in my_A_p
  bordered_grp->extractParameterComponent(false, *my_A, my_A_p);

  // Extract solution component from my_B and store in my_B_p
  if (isZeroB)
    my_B_p.putScalar(0.0);
  else
    bordered_grp->extractParameterComponent(true, *my_B, my_B_p);

  // Copy in my_C
  my_CC.assign(*my_C);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeInverse(
                      Teuchos::ParameterList& params,
                      const NOX::Abstract::Vector& input,
                      NOX::Abstract::Vector& result) const
{
  // Convert input, result to multivectors
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_input =
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_result =
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianTransposeInverse
  NOX::Abstract::Group::ReturnType status =
    applyJacobianTransposeInverseMultiVector(params, *mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeInverseMultiVector(
                     Teuchos::ParameterList& params,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeInverseMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

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

  // Call bordered solver applyInverseTranspose method
  status =
    borderedSolver->initForTransposeSolve();
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                finalStatus,
                                callingFunction);
  status = borderedSolver->applyInverseTranspose(params, input_x.get(),
                         input_param.get(),
                         *result_x, *result_param);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                finalStatus,
                                callingFunction);

  return finalStatus;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ConstrainedGroup::getGroup()
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::ConstrainedGroup::getConstraints()
{
  return constraintsPtr;
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
