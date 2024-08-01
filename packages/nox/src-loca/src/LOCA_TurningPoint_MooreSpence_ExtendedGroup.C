// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_TurningPoint_MooreSpence_ExtendedGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_TurningPoint_MooreSpence_SolverStrategy.H"
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "NOX_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_TimeDependent_AbstractGroup.H"

LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
     const Teuchos::RCP<LOCA::TurningPoint::MooreSpence::AbstractGroup>& g)
  : LOCA::Extended::MultiAbstractGroup(),
    LOCA::MultiContinuation::AbstractGroup(),
    globalData(global_data),
    parsedParams(topParams),
    turningPointParams(tpParams),
    grpPtr(g),
    xMultiVec(globalData, g->getX(), 1),
    fMultiVec(globalData, g->getX(), 2),
    newtonMultiVec(globalData, g->getX(), 1),
    lengthMultiVec(),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    lengthVec(),
    solverStrategy(),
    index_f(1),
    index_dfdp(1),
    bifParamID(1),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false),
    updateVectorsEveryContinuationStep(false),
    nullVecScaling(NVS_OrderN),
    multiplyMass(false),
    tdGrp(),
    tmp_mass()
{
  const char *func = "LOCA::TurningPoint::MooreSpence::ExtendedGroup()";

  // Set x
  *(xMultiVec.getColumn(0)->getXVec()) = g->getX();

  if (!turningPointParams->isParameter("Bifurcation Parameter")) {
    globalData->locaErrorCheck->throwError(func,
                 "\"Bifurcation Parameter\" name is not set!");
  }
  std::string bifParamName = turningPointParams->get(
                          "Bifurcation Parameter",
                          "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamID[0] = p.getIndex(bifParamName);

  if (!turningPointParams->isParameter("Length Normalization Vector")) {
    globalData->locaErrorCheck->throwError(func,
               "\"Length Normalization Vector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> lenVecPtr =
    (*turningPointParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<NOX::Abstract::Vector> >("Length Normalization Vector");

  if (!turningPointParams->isParameter("Initial Null Vector")) {
    globalData->locaErrorCheck->throwError(func,
                 "\"Initial Null Vector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> nullVecPtr =
    (*turningPointParams).INVALID_TEMPLATE_QUALIFIER
    get<Teuchos::RCP<NOX::Abstract::Vector> >("Initial Null Vector");

  bool perturbSoln = turningPointParams->get(
                           "Perturb Initial Solution",
                           false);
  double perturbSize = turningPointParams->get(
                         "Relative Perturbation Size",
                         1.0e-3);

  updateVectorsEveryContinuationStep =
    turningPointParams->get("Update Null Vectors Every Continuation Step",
                false);
  std::string nullVecScalingMethod =
    turningPointParams->get("Null Vector Scaling", "Order N");
  if (nullVecScalingMethod == "None")
    nullVecScaling = NVS_None;
  else if (nullVecScalingMethod == "Order 1")
    nullVecScaling = NVS_OrderOne;
  else if (nullVecScalingMethod == "Order N")
    nullVecScaling = NVS_OrderN;
  else
    globalData->locaErrorCheck->throwError(
       "LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup()",
       std::string("Unknown null vector scaling method:  ") + nullVecScalingMethod);
   multiplyMass =
    turningPointParams->get("Multiply Null Vectors by Mass Matrix", false);
  if (multiplyMass && tdGrp == Teuchos::null) {
    globalData->locaErrorCheck->throwError(
       "LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup()",
       "Group must be derived from LOCA::TimeDependent::AbstractGroup to multiply null vectors by mass matrix");
  }

  lengthMultiVec =
    lenVecPtr->createMultiVector(1, NOX::DeepCopy);
  *(xMultiVec.getColumn(0)->getNullVec()) = *nullVecPtr;

  tdGrp = Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grpPtr);
  tmp_mass = grpPtr->getX().clone(NOX::ShapeCopy);

  // Instantiate solver strategy
  solverStrategy =
    globalData->locaFactory->createMooreSpenceTurningPointSolverStrategy(
                              parsedParams,
                              turningPointParams);

  // Set up multi-vector views
  setupViews();

  init(perturbSoln, perturbSize);
}

LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup(
        const LOCA::TurningPoint::MooreSpence::ExtendedGroup& source,
        NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    turningPointParams(source.turningPointParams),
    grpPtr(Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MooreSpence::AbstractGroup>(source.grpPtr->clone(type))),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    lengthMultiVec(source.lengthMultiVec->clone(type)),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    lengthVec(),
    solverStrategy(source.solverStrategy),
    index_f(1),
    index_dfdp(1),
    bifParamID(source.bifParamID),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton),
    updateVectorsEveryContinuationStep(source.updateVectorsEveryContinuationStep),
    nullVecScaling(source.nullVecScaling),
    multiplyMass(source.multiplyMass),
    tdGrp(source.tdGrp),
    tmp_mass(source.tmp_mass->clone(type))
{

  // Instantiate solver strategy
  solverStrategy =
    globalData->locaFactory->createMooreSpenceTurningPointSolverStrategy(
                                   parsedParams,
                                   turningPointParams);

  // Set up multi-vector views
  setupViews();

  if (type == NOX::ShapeCopy) {
    isValidF = false;
    isValidJacobian = false;
    isValidNewton = false;
  }
}

LOCA::TurningPoint::MooreSpence::ExtendedGroup::~ExtendedGroup()
{
}

NOX::Abstract::Group&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::operator=(
                       const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::TurningPoint::MooreSpence::ExtendedGroup::clone(
                            NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedGroup(*this,
                                    type));
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setX(
                          const NOX::Abstract::Vector& y)
{
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& yy =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(y);
  grpPtr->setX( *yy.getXVec() );
  *xVec = y;
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeX(
                          const NOX::Abstract::Group& g,
                          const NOX::Abstract::Vector& d,
                          double step)
{
  const LOCA::TurningPoint::MooreSpence::ExtendedGroup& gg =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedGroup&>(g);
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& dd =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(d);

  grpPtr->computeX(*(gg.grpPtr), *dd.getXVec(), step);
  xVec->update(1.0, gg.getX(), step, dd, 0.0);
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeF()
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                           callingFunction);
  }
  *(fVec->getXVec()) = grpPtr->getF();

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                           callingFunction);
  }

  // Compute J*n
  status = grpPtr->applyJacobian(*(xVec->getNullVec()),
                 *(fVec->getNullVec()));
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  // Compute phi^T*n
  fVec->getBifParam() = lTransNorm(*(xVec->getNullVec())) - 1.0;

  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeJacobian()
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  status = grpPtr->computeDfDpMulti(bifParamID,
                    *fMultiVec.getXMultiVec(),
                    isValidF);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  // Compute underlying dJn/dp (may invalidate underlying data)
  status = grpPtr->computeDJnDpMulti(bifParamID,
                     *(xVec->getNullVec()),
                     *fMultiVec.getNullMultiVec(),
                     isValidF);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  // Compute underlying Jacobian
  status = grpPtr->computeJacobian();
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  solverStrategy->setBlocks(
          grpPtr,
          Teuchos::rcp(this, false),
          xVec->getNullVec(),
          fVec->getNullVec(),
          dfdpMultiVec->getXMultiVec(),
          dfdpMultiVec->getNullMultiVec());

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeGradient()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeNewton(
                         Teuchos::ParameterList& params)
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeNewton()";
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

  // solve using contiguous
  status = solverStrategy->solve(params, *ffMultiVec, newtonMultiVec);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  newtonMultiVec.scale(-1.0);

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobian(
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
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTranspose(
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
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverse(
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
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianMultiVector(
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
                 "Called with invalid Jacobian!");
  }

  // Cast vectors to TP vectors
  const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_input =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_result =
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(result);

  // Get constant references to input vector components
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x =
    tp_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_null =
    tp_input.getNullMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param =
    tp_input.getScalars();

  // Get non-constant references to result vector components
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x =
    tp_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_null =
    tp_result.getNullMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param =
    tp_result.getScalars();

  // Temporary vector
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp =
    input_null->clone(NOX::ShapeCopy);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // compute J*x
  status = grpPtr->applyJacobianMultiVector(*input_x, *result_x);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // compute J*x + p*dR/dp
  result_x->update(Teuchos::NO_TRANS, 1.0, *(dfdpMultiVec->getXMultiVec()),
           *input_param, 1.0);

  // compute J*y
  status = grpPtr->applyJacobianMultiVector(*input_null, *result_null);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // compute J*y + p*dJy/dp
  result_null->update(Teuchos::NO_TRANS, 1.0,
              *(dfdpMultiVec->getNullMultiVec()),
              *input_param, 1.0);

  // compute (dJy/dx)*x
  status = grpPtr->computeDJnDxaMulti(*(xVec->getNullVec()),
                      *(fVec->getNullVec()),
                      *input_x, *tmp);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null->update(1.0, *tmp, 1.0);

  // compute l^T*y
  lTransNorm(*input_null, *result_param);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector(
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
                 "Called with invalid Jacobian!");
  }

  // Cast vectors to TP vectors
  const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_input =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_result =
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(result);

  // Get constant references to input vector components
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x =
    tp_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_null =
    tp_input.getNullMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param =
    tp_input.getScalars();

  // Get non-constant references to result vector components
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x =
    tp_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_null =
    tp_result.getNullMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param =
    tp_result.getScalars();

  // Temporary vector
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp =
    input_null->clone(NOX::ShapeCopy);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // compute J^T*x
  status = grpPtr->applyJacobianTransposeMultiVector(*input_x, *result_x);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // compute J^T*y
  status = grpPtr->applyJacobianTransposeMultiVector(*input_null, *result_null);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // compute (dJn/dx)^T*y
  status = grpPtr->computeDwtJnDxMulti(*input_null, *(xVec->getNullVec()),
                       *tmp);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // compute J^T*x + (dJn/dx)^T*y
  result_x->update(1.0, *tmp, 1.0);

  // compute J^T*y + phi*z
  result_null->update(Teuchos::NO_TRANS, 1.0/lengthVec->length(), *lengthMultiVec, *input_param, 1.0);

  // compute (df/dp)^T*x
  input_x->multiply(1.0, *(dfdpMultiVec->getXMultiVec()), *result_param);

  // compute (dJn/dp)^T*y
  NOX::Abstract::MultiVector::DenseMatrix t(1, input_param->numCols());
  input_null->multiply(1.0, *(dfdpMultiVec->getNullMultiVec()), t);

  // compute (df/dp)^T*x + (dJn/dp)^T*y
  *result_param += t;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector(
                     Teuchos::ParameterList& params,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  // Cast vectors to TP vectors
  const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_input =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_result =
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(result);

  return solverStrategy->solve(params, tp_input, tp_result);
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isF() const
{
  return isValidF;
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isJacobian() const
{
  return isValidJacobian;
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isGradient() const
{
  return false;
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isNewton() const
{
  return isValidNewton;
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getX() const
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getF() const
{
  return *fVec;
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormF() const
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradient() const
{
  globalData->locaErrorCheck->throwError(
          "LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradient()",
          " - not implemented");
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNewton() const
{
  return *newtonVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getXPtr() const
{
  return xVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getFPtr() const
{
  return fVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradientPtr() const
{
  globalData->locaErrorCheck->throwError(
          "LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradientPtr()",
          " - not implemented");
  return getNewtonPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNewtonPtr() const
{
  return newtonVec;
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual() const
{
  std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::TurningPoint::MooreSpence::ExtendedVector residual = *fVec;

  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup()
{
  return grpPtr;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::copy(
                        const NOX::Abstract::Group& src)
{
  const LOCA::TurningPoint::MooreSpence::ExtendedGroup& source =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedGroup&>(src);

  // Protect against A = A
  if (this != &source) {

    // Copy values
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    turningPointParams = source.turningPointParams;
    grpPtr->copy(*(source.grpPtr));
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    *lengthMultiVec = *source.lengthMultiVec;
    index_f = source.index_f;
    index_dfdp = source.index_dfdp;
    bifParamID = source.bifParamID;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    updateVectorsEveryContinuationStep =
      source.updateVectorsEveryContinuationStep;
    nullVecScaling = source.nullVecScaling;
    multiplyMass = source.multiplyMass;

    // set up views again just to be safe
    setupViews();

    // Instantiate solver strategy
    solverStrategy =
      globalData->locaFactory->createMooreSpenceTurningPointSolverStrategy(
                                           parsedParams,
                                           turningPointParams);
  }
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParamsMulti(
              const std::vector<int>& paramIDs,
              const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  for (unsigned int i=0; i<paramIDs.size(); i++) {
    if (paramIDs[i] == bifParamID[0])
      setBifParam(vals(i,0));
  }
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParams(
                          const LOCA::ParameterVector& p)
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
  setBifParam(p[bifParamID[0]]);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam(int paramID,
                             double val)
{
  if (paramID == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam(std::string paramID,
                             double val)
{
  const LOCA::ParameterVector& pVec = grpPtr->getParams();
  if (pVec.getIndex(paramID) == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

const LOCA::ParameterVector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParams() const
{
  return grpPtr->getParams();
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam(std::string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDpMulti(
                        const std::vector<int>& paramIDs,
                        NOX::Abstract::MultiVector& dfdp,
                        bool isValid_F)
{
   std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDpMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast dfdp to TP vector
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_dfdp =
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(dfdp);

  // Compute df/dp
  status = grpPtr->computeDfDpMulti(paramIDs, *tp_dfdp.getXMultiVec(),
                    isValid_F);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Compute d(Jn)/dp
  status = grpPtr->computeDJnDpMulti(paramIDs, *(xVec->getNullVec()),
                     *tp_dfdp.getNullMultiVec(), isValid_F);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Set parameter components
  if (!isValid_F)
    tp_dfdp.getScalar(0,0) = lTransNorm(*(xVec->getNullVec()));
  for (int i=0; i<dfdp.numVectors()-1; i++)
    tp_dfdp.getScalar(0,i+1) = 0.0;

  return finalStatus;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::preProcessContinuationStep(
             LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  // Rescale length vector so that the normalization condition is met
  double lVecDotNullVec = lTransNorm(*(xVec->getNullVec()));

  if (lVecDotNullVec == 0.0) {
    globalData->locaErrorCheck->throwError(
      "LOCA::TurningPoint::MooreSpence::ExtendedGroup::preProcessContinuationStep()",
      "null vector can be orthogonal to length-scaling vector");
  }
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
    "\tIn LOCA::TurningPoint::MooreSpence::ExtendedGroup::preProcessContinuationStep(), " <<
    "scaling null vector by:" <<
    globalData->locaUtils->sciformat(1.0 / lVecDotNullVec) << std::endl;
  }
  if (updateVectorsEveryContinuationStep) {
    xVec->getNullVec()->scale(1.0/lVecDotNullVec);
  }

  grpPtr->preProcessContinuationStep(stepStatus);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::postProcessContinuationStep(
             LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful &&
      updateVectorsEveryContinuationStep) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
      globalData->locaUtils->out() <<
      "\n\tUpdating null vector for the next continuation step" << std::endl;
    }
    *lengthVec = *(xVec->getNullVec());

    scaleNullVector(*lengthVec);
  }

  grpPtr->postProcessContinuationStep(stepStatus);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDraw(
                           const NOX::Abstract::Vector& x,
                           double *px) const
{
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& mx =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(x);

  grpPtr->projectToDraw(*(mx.getXVec()), px);
  px[grpPtr->projectToDrawDimension()] = mx.getBifParam();
}

int
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension() + 1;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution(
                          const double conParam) const
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Turning Point located at: " <<
      globalData->locaUtils->sciformat(conParam) << "   " <<
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;

    globalData->locaUtils->out() <<
      "\tPrinting Solution Vector for conParam = " <<
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tPrinting Null Vector for bif param = " <<
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*(xVec->getNullVec()), xVec->getBifParam());
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution(
                         const NOX::Abstract::Vector& x_,
                         const double conParam) const
{
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& tp_x =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(x_);

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Turning Point located at: " <<
      globalData->locaUtils->sciformat(conParam) << "   " <<
      globalData->locaUtils->sciformat(tp_x.getBifParam()) << std::endl;

     globalData->locaUtils->out() <<
       "\tPrinting Solution Vector for conParam = " <<
       globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(*tp_x.getXVec(), conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tPrinting Null Vector for bif param = " <<
      globalData->locaUtils->sciformat(tp_x.getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*tp_x.getNullVec(), tp_x.getBifParam());
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeInverse(
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
    applyJacobianTransposeInverseMultiVector(params, *mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeInverseMultiVector(
                     Teuchos::ParameterList& params,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  // Cast vectors to TP vectors
  const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_input =
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_result =
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(result);

  return solverStrategy->solveTranspose(params, tp_input, tp_result);
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getBifParam() const
{
  return grpPtr->getParam(bifParamID[0]);
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm(
                    const NOX::Abstract::Vector& n) const
{
  return lengthVec->innerProduct(n) / lengthVec->length();
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm(
            const NOX::Abstract::MultiVector& n,
            NOX::Abstract::MultiVector::DenseMatrix& result) const
{
  n.multiply(1.0 / lengthVec->length(), *lengthMultiVec, result);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getLengthVector() const
{
  Teuchos::RCP<NOX::Abstract::Vector> l = lengthVec->clone();
  l->scale(1.0 / lengthVec->length());
  return l;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setBifParam(double param)
{
  grpPtr->setParam(bifParamID[0], param);
  xVec->getBifParam() = param;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setupViews()
{
  index_f[0] = 0;
  index_dfdp[0] = 1;

  xVec = xMultiVec.getColumn(0);
  fVec = fMultiVec.getColumn(0);
  newtonVec = newtonMultiVec.getColumn(0);
  lengthVec = Teuchos::rcp(&(*lengthMultiVec)[0],false);

  ffMultiVec = Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector>(fMultiVec.subView(index_f),true);

  dfdpMultiVec = Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector>(fMultiVec.subView(index_dfdp),true);

}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::init(bool perturbSoln,
                             double perturbSize)
{
  xVec->getBifParam() = getBifParam();

  // Rescale null vector
  scaleNullVector(*lengthVec);

  // Rescale length vector so that the normalization condition is met
  double lVecDotNullVec = lTransNorm(*(xVec->getNullVec()));

  if (lVecDotNullVec == 0.0) {
    globalData->locaErrorCheck->throwError(
           "LOCA::TurningPoint::MooreSpence::ExtendedGroup::init()",
           "null vector can be orthogonal to length-scaling vector");
  }
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tIn LOCA::TurningPoint::MooreSpence::ExtendedGroup::init(), " <<
      "scaling null vector by:" <<
      globalData->locaUtils->sciformat(1.0 / lVecDotNullVec) << std::endl;
  }
  xVec->getNullVec()->scale(1.0/lVecDotNullVec);

  if (perturbSoln) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
     globalData->locaUtils->out() <<
       "\tIn LOCA::TurningPoint::MooreSpence::ExtendedGroup::init(), " <<
       "applying random perturbation to initial solution of size: " <<
       globalData->locaUtils->sciformat(perturbSize) << std::endl;
    }
    Teuchos::RCP<NOX::Abstract::Vector> perturb =
      xVec->getXVec()->clone(NOX::ShapeCopy);
    perturb->random();
    perturb->scale(*(xVec->getXVec()));
    xVec->getXVec()->update(perturbSize, *perturb, 1.0);
    grpPtr->setX(*(xVec->getXVec()));
  }
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::
scaleNullVector(NOX::Abstract::Vector& a)
{
  std::string callingFunction =
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::scaleNullVector()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  if (multiplyMass && tdGrp != Teuchos::null) {
    status = tdGrp->computeShiftedMatrix(0.0, 1.0);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
    *tmp_mass = a;
    status = tdGrp->applyShiftedMatrix(*tmp_mass, a);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
    status = tdGrp->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }
  if (nullVecScaling == NVS_OrderOne) {
    a.scale(1.0 / a.norm());
  }
  else if (nullVecScaling == NVS_OrderN) {
    double dn = a.length();
    a.scale(std::sqrt(dn) / a.norm());
  }
}
