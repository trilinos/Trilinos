// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Pitchfork_MinimallyAugmented_ExtendedGroup.H"

#include "Teuchos_ParameterList.hpp"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_Constraint.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"
#include "LOCA_BorderedSolver_JacobianOperator.H"

LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
ExtendedGroup(
       const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& pfParams,
       const Teuchos::RCP<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>& g)
  : globalData(global_data),
    parsedParams(topParams),
    pitchforkParams(pfParams),
    grpPtr(g),
    bordered_grp(),
    constraintsPtr(),
    xMultiVec(globalData, g->getX(), 1, 2, NOX::DeepCopy),
    fMultiVec(globalData, g->getX(), 3, 2, NOX::ShapeCopy),
    newtonMultiVec(globalData, g->getX(), 1, 2, NOX::ShapeCopy),
    gradientMultiVec(globalData, g->getX(), 1, 2, NOX::ShapeCopy),
    xVec(),
    psiVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    fBifMultiVec(),
    newtonVec(),
    gradientVec(),
    jacOp(),
    borderedSolver(),
    index_f(1),
    index_dfdp(2),
    bifParamID(-1),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false),
    isValidGradient(false),
    isBordered(false)
{
  const char *func = "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup()";

  // Set up multi-vector views
  setupViews();

  // Get bifurcation parameter name
  if (!pitchforkParams->isParameter("Bifurcation Parameter")) {
    globalData->locaErrorCheck->throwError(func,
                 "\"Bifurcation Parameter\" name is not set!");
  }
  std::string bifParamName = pitchforkParams->get("Bifurcation Parameter",
                         "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamID = p.getIndex(bifParamName);

  // Get psi vector
  if (!pitchforkParams->isParameter("Antisymmetric Vector")) {
    globalData->locaErrorCheck->throwError(func,
               "\"Antisymmetric Vector\" is not set!");
  }
  psiVec =
    (*pitchforkParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<NOX::Abstract::Vector> >("Antisymmetric Vector");

  // Create constraint equation
  constraintsPtr =
    Teuchos::rcp(new LOCA::Pitchfork::MinimallyAugmented::Constraint(
                                   globalData,
                                   parsedParams,
                                   pfParams,
                                   grpPtr,
                                   psiVec,
                                   bifParamID));

  // Set parameters in solution vector
  xVec->getScalar(0) = grpPtr->getParam(bifParamID);
  xVec->getScalar(1) = 0.0;   // initial value of slack variable

  // Set parameters and solution vector in constraints
  constraintsPtr->setParam(bifParamID, xVec->getScalar(0));
  constraintsPtr->setX(*(xVec->getXVec()));

  // Instantiate bordered solver
  borderedSolver = globalData->locaFactory->createBorderedSolverStrategy(
                   parsedParams,
                   pitchforkParams);

  // Determine if underlying group is bordered
  bordered_grp =
    Teuchos::rcp_dynamic_cast<LOCA::BorderedSystem::AbstractGroup>(grpPtr);
  isBordered = (bordered_grp != Teuchos::null);

  // Create Jacobian operator for bordered solver
  jacOp = Teuchos::rcp(new LOCA::BorderedSolver::JacobianOperator(grpPtr));
}

LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
ExtendedGroup(const LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup& source,
          NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    pitchforkParams(source.pitchforkParams),
    grpPtr(Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>(source.grpPtr->clone(type))),
    bordered_grp(),
    constraintsPtr(Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MinimallyAugmented::Constraint>(source.constraintsPtr->clone(type))),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    gradientMultiVec(source.gradientMultiVec, type),
    xVec(),
    psiVec(source.psiVec),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    fBifMultiVec(),
    newtonVec(),
    gradientVec(),
    jacOp(),
    borderedSolver(source.borderedSolver),
    index_f(1),
    index_dfdp(2),
    bifParamID(source.bifParamID),
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
                   pitchforkParams);

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

  constraintsPtr->setGroup(grpPtr);

  // Create Jacobian operator for bordered solver
  jacOp = Teuchos::rcp(new LOCA::BorderedSolver::JacobianOperator(grpPtr));

  // Set blocks in bordered solver
  if (isValidJacobian) {
    borderedSolver->setMatrixBlocks(jacOp,
                    dfdpMultiVec->getXMultiVec(),
                    constraintsPtr,
                    dfdpMultiVec->getScalars());
    NOX::Abstract::Group::ReturnType status = borderedSolver->initForSolve();
    globalData->locaErrorCheck->checkReturnType(status,
                        "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup()");
  }
}


LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
~ExtendedGroup()
{
}

double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getBifParam() const
{
  return grpPtr->getParam(bifParamID);
}

NOX::Abstract::Group&
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
operator=(const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new ExtendedGroup(*this, type));
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
setX(const NOX::Abstract::Vector& y)
{
  const LOCA::MultiContinuation::ExtendedVector& my =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  grpPtr->setX( *(my.getXVec()) );
  grpPtr->setParam(bifParamID, my.getScalar(0));
  *xVec = my;
  constraintsPtr->setX( *(my.getXVec()) );
  constraintsPtr->setParam(bifParamID, my.getScalar(0));

  resetIsValid();
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
computeX(const NOX::Abstract::Group& g,
     const NOX::Abstract::Vector& d,
     double step)
{
  const LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup& mg =
    dynamic_cast<const LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup&>(g);
  const LOCA::MultiContinuation::ExtendedVector& md =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(d);

  grpPtr->computeX(*(mg.grpPtr), *(md.getXVec()), step);
  xVec->update(1.0, mg.getX(), step, md, 0.0);
  grpPtr->setParam(bifParamID, xVec->getScalar(0));
  constraintsPtr->setX( *(xVec->getXVec()) );
  constraintsPtr->setParam(bifParamID, xVec->getScalar(0));

  resetIsValid();
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
computeF()
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeF()";
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
  fVec->getXVec()->update(1.0, grpPtr->getF(), xVec->getScalar(1), *psiVec,
              0.0);

  // Compute constraints
  if (!constraintsPtr->isConstraints()) {
    status = constraintsPtr->computeConstraints();
  }
  fVec->getScalars()->assign(constraintsPtr->getConstraints());

  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
computeJacobian()
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  // We force recomputation of f since we store a different f than the
  // underlying group
  std::vector<int> paramIDs(1);
  paramIDs[0] = bifParamID;
  status = grpPtr->computeDfDpMulti(paramIDs,
                    *(fBifMultiVec->getXMultiVec()),
                    false);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Add on slack component
  fVec->getXVec()->update(xVec->getScalar(1), *psiVec, 1.0);

  // Compute derivative w.r.t. slack variable
  (*(dfdpMultiVec->getXMultiVec()))[1] = *psiVec;

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
  status = constraintsPtr->computeDP(paramIDs,
                     *(fBifMultiVec->getScalars()),
                     isValidF);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                finalStatus,
                                callingFunction);

  // Compute derivative w.r.t. slack variable
  dfdpMultiVec->getScalar(0,1) = 0.0;
  dfdpMultiVec->getScalar(1,1) = 0.0;

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Set blocks in bordered solver
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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
computeGradient()
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeGradient()";
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

  // Compute J^T*f for pitchfork group
  status = applyJacobianTranspose(*fVec, *gradientVec);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  isValidGradient = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
computeNewton(Teuchos::ParameterList& params)
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeNewton()";
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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
applyJacobianMultiVector(const NOX::Abstract::MultiVector& input,
             NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianMultiVector()";

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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
applyJacobianTransposeMultiVector(const NOX::Abstract::MultiVector& input,
                  NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeMultiVector()";

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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
applyJacobianInverseMultiVector(Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& input,
                NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianInverseMultiVector()";

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
  NOX::Abstract::Group::ReturnType status =
    borderedSolver->applyInverse(params, input_x.get(), input_param.get(),
                 *result_x, *result_param);

  return status;
}

bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
isF() const
{
  return isValidF;
}

bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
isJacobian() const
{
  return isValidJacobian;
}

bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
isGradient() const
{
  return isValidGradient;
}

bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
isNewton() const
{
  return isValidNewton;
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getX() const
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getF() const
{
  return *fVec;
}

double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getNormF() const
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getGradient() const
{
  return *gradientVec;
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getNewton() const
{
  return *newtonVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getXPtr() const
{
  return xVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getFPtr() const
{
  return fVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getGradientPtr() const
{
  return gradientVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getNewtonPtr() const
{
  return newtonVec;
}

double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getNormNewtonSolveResidual() const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::MultiContinuation::ExtendedVector residual = *fVec;

  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual = residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getUnderlyingGroup()
{
  return grpPtr;
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
copy(const NOX::Abstract::Group& src)
{

  const LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup& source =
    dynamic_cast<const LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup&>(src);

  // Protect against A = A
  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    pitchforkParams = source.pitchforkParams;
    grpPtr->copy(*source.grpPtr);
    constraintsPtr->copy(*source.constraintsPtr);
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    gradientMultiVec = source.gradientMultiVec;
    index_f = source.index_f;
    index_dfdp = source.index_dfdp;
    bifParamID = source.bifParamID;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;

    // set up views again just to be safe
    setupViews();

    // Instantiate bordered solver
    borderedSolver =
      globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
                                pitchforkParams);

    // Set blocks in bordered solver
    if (isValidJacobian) {
      borderedSolver->setMatrixBlocks(jacOp,
                      dfdpMultiVec->getXMultiVec(),
                      constraintsPtr,
                      dfdpMultiVec->getScalars());
      NOX::Abstract::Group::ReturnType status = borderedSolver->initForSolve();
      globalData->locaErrorCheck->checkReturnType(status,
                          "LOCA::Pitchfork::MinimallyAugmented::copy()");
  }
  }
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
setParamsMulti(const std::vector<int>& paramIDs,
           const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  constraintsPtr->setParams(paramIDs, vals);

  for (unsigned int i=0; i<paramIDs.size(); i++)
    if (paramIDs[i] == bifParamID)
    xVec->getScalar(0) = vals(i,0);

  resetIsValid();
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
setParams(const LOCA::ParameterVector& p)
{
  grpPtr->setParams(p);
  for (int i=0; i<p.length(); i++)
    constraintsPtr->setParam(i, p[i]);
  xVec->getScalar(0) = p[bifParamID];

  resetIsValid();
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
  constraintsPtr->setParam(paramID, val);

  if (paramID == bifParamID)
    xVec->getScalar(0) = val;

  resetIsValid();
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
setParam(std::string paramID, double val)
{
  const LOCA::ParameterVector& p = grpPtr->getParams();
  int id = p.getIndex(paramID);
  setParam(id, val);
}

const LOCA::ParameterVector&
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getParams() const
{
  return grpPtr->getParams();
}

double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getParam(std::string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
computeDfDpMulti(const std::vector<int>& paramIDs,
         NOX::Abstract::MultiVector& dfdp,
         bool isValid_F)
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeDfDpMulti()";
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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->preProcessContinuationStep(stepStatus);
  constraintsPtr->preProcessContinuationStep(stepStatus);
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->postProcessContinuationStep(stepStatus);
  constraintsPtr->postProcessContinuationStep(stepStatus);
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
projectToDraw(const NOX::Abstract::Vector& x,
          double *px) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);

  grpPtr->projectToDraw(*mx.getXVec(), px);
  for (int i=0; i<2; i++)
    px[grpPtr->projectToDrawDimension()+i] = mx.getScalar(i);
}

int
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension() + 2;
}

double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
computeScaledDotProduct(const NOX::Abstract::Vector& a,
            const NOX::Abstract::Vector& b) const
{
  const LOCA::MultiContinuation::ExtendedVector& ma =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(a);
  const LOCA::MultiContinuation::ExtendedVector& mb =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(b);

  double val = grpPtr->computeScaledDotProduct(*ma.getXVec(), *mb.getXVec());
  for (int i=0; i<2; i++) {
    val += ma.getScalar(i) * mb.getScalar(i);
  }

  return val;
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
printSolution(const double conParam) const
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Pitchfork located at: " <<
      globalData->locaUtils->sciformat(conParam) << "   " <<
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;

    globalData->locaUtils->out() <<
      "\tSlack variable = " <<
      globalData->locaUtils->sciformat(xVec->getScalar(1)) << std::endl;

    globalData->locaUtils->out() <<
      "\tPrinting Solution Vector for conParam = " <<
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tPrinting Right Null Vector for bif param = " <<
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*(constraintsPtr->getRightNullVec()), getBifParam());
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tPrinting Left Null Vector for sigma = " <<
      globalData->locaUtils->sciformat(constraintsPtr->getSigma()) <<
      std::endl;
  }
  grpPtr->printSolution(*(constraintsPtr->getLeftNullVec()),
            constraintsPtr->getSigma());
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
printSolution(const NOX::Abstract::Vector& x,
          const double conParam) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Pitchfork located at: " <<
      globalData->locaUtils->sciformat(conParam) << "   " <<
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;

    globalData->locaUtils->out() <<
      "\tSlack variable = " <<
      globalData->locaUtils->sciformat(mx.getScalar(1)) << std::endl;

    globalData->locaUtils->out() <<
      "\tPrinting Solution Vector for conParam = " <<
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(*mx.getXVec(), conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tPrinting Right Null Vector for bif param = " <<
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*(constraintsPtr->getRightNullVec()), getBifParam());
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tPrinting Left Null Vector for sigma = " <<
      globalData->locaUtils->sciformat(constraintsPtr->getSigma()) <<
      std::endl;
  }
  grpPtr->printSolution(*(constraintsPtr->getLeftNullVec()),
            constraintsPtr->getSigma());
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
scaleVector(NOX::Abstract::Vector& x) const
{
  LOCA::MultiContinuation::ExtendedVector& mx =
    dynamic_cast<LOCA::MultiContinuation::ExtendedVector&>(x);

  grpPtr->scaleVector(*mx.getXVec());
}

int
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getBorderedWidth() const
{
  int my_width = 2;
  if (isBordered)
    return my_width + bordered_grp->getBorderedWidth();
  else
    return my_width;
}

Teuchos::RCP<const NOX::Abstract::Group>
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
getUnborderedGroup() const
{
  if (isBordered)
    return bordered_grp->getUnborderedGroup();
  else
    return grpPtr;
}

bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
isCombinedAZero() const
{
  return false;  // A is always considered non-zero (A = [df/dp psi])
}

bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
isCombinedBZero() const
{
  return false; // B is always considered non-zero (B = [dsigma/dx psi])
}

bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
isCombinedCZero() const
{
  return false;  // C is always considered non-zero (C = [dsigma/dp 0; 0 0] )
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
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
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
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
                          2, num_cols, w, 0);

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
                          num_rows, 2, 0, w);

    // Decompose mc_v_x
    bordered_grp->extractParameterComponent(use_transpose,*mc_v_x, v_p_1);
    for (int j=0; j<2; j++)
      for (int i=0; i<num_rows; i++)
    v_p_2(i,j) = (*mc_v_p)(j,i);
  }
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
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
                        2, num_cols, w, 0);

  // load v_x, v_p_1 into mc_v_x
  bordered_grp->loadNestedComponents(v_x, v_p_1, *mc_v_x);

  // load v_p_2 into mc_v_p
  mc_v_p->assign(v_p_2);
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
fillA(NOX::Abstract::MultiVector& A) const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillA";

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
  std::vector<int> idx2(2);
  for (int i=0; i<2; i++)
    idx2[i] = w+i;
  Teuchos::RCP<NOX::Abstract::MultiVector> my_A_x =
    A.subView(idx2);

  // Extract solution component from my_A and store in A
  bordered_grp->extractSolutionComponent(*my_A, *my_A_x);
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
fillB(NOX::Abstract::MultiVector& B) const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillB";

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_B =
    Teuchos::rcp(constraintsPtr->getDX(),false);

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
  for (int i=0; i<2; i++)
    idx2[i] = w+i;
  Teuchos::RCP<NOX::Abstract::MultiVector> my_B_x =
    B.subView(idx2);

  // Extract solution component from my_B and store in B
  bordered_grp->extractSolutionComponent(*my_B, *my_B_x);
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
fillC(NOX::Abstract::MultiVector::DenseMatrix& C) const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillC";

  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> my_C =
    dfdpMultiVec->getScalars();

  // If the underlying system isn't bordered, we're done
  if (!isBordered) {
    C.assign(*my_C);
    return;
  }

  Teuchos::RCP<const NOX::Abstract::MultiVector> my_B =
    Teuchos::rcp(constraintsPtr->getDX(),false);

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
                         w, 2, 0, w);
  NOX::Abstract::MultiVector::DenseMatrix my_B_p(Teuchos::View, C,
                         2, w, w, 0);
  NOX::Abstract::MultiVector::DenseMatrix my_CC(Teuchos::View, C,
                        2, 2, w, w);

  // Extract solution component from my_A and store in my_A_p
  bordered_grp->extractParameterComponent(false, *my_A, my_A_p);

  // Extract solution component from my_B and store in my_B_p
  bordered_grp->extractParameterComponent(true, *my_B, my_B_p);

  // Copy in my_C
  my_CC.assign(*my_C);
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
resetIsValid()
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
setupViews()
{
  index_f[0] = 0;
  for (int i=0; i<2; i++)
    index_dfdp[i] = i+1;

  xVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(xMultiVec.getVector(0),true);
  fVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(fMultiVec.getVector(0),true);
  newtonVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(newtonMultiVec.getVector(0),true);
  gradientVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(gradientMultiVec.getVector(0),true);

  ffMultiVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(fMultiVec.subView(index_f),true);

  dfdpMultiVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(fMultiVec.subView(index_dfdp),true);

  std::vector<int> index_fbif(2);
  index_fbif[0] = 0; index_fbif[1] = 1;
  fBifMultiVec = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(fMultiVec.subView(index_fbif),true);

}

void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::
setBifParam(double val)
{
  grpPtr->setParam(bifParamID, val);
  xVec->getScalar(0) = val;
  constraintsPtr->setParam(bifParamID, val);

  resetIsValid();
}

