// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Pitchfork_MooreSpence_ExtendedGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_Pitchfork_MooreSpence_SolverStrategy.H"
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "NOX_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::Pitchfork::MooreSpence::ExtendedGroup::ExtendedGroup(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
     const Teuchos::RCP<LOCA::Pitchfork::MooreSpence::AbstractGroup>& g)
  : LOCA::Extended::MultiAbstractGroup(),
    LOCA::MultiContinuation::AbstractGroup(),
    globalData(global_data),
    parsedParams(topParams),
    pitchforkParams(tpParams),
    grpPtr(g),
    xMultiVec(globalData, g->getX(), 1),
    fMultiVec(globalData, g->getX(), 2),
    newtonMultiVec(globalData, g->getX(), 1),
    asymMultiVec(),
    lengthMultiVec(),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    asymVec(),
    lengthVec(),
    solverStrategy(),
    index_f(1),
    index_dfdp(1),
    bifParamID(1),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  const char *func = "LOCA::Pitchfork::MooreSpence::ExtendedGroup()";

  // Set x
  *(xMultiVec.getColumn(0)->getXVec()) = g->getX();

  if (!pitchforkParams->isParameter("Bifurcation Parameter")) {
    globalData->locaErrorCheck->throwError(func,
                 "\"Bifurcation Parameter\" name is not set!");
  }
  std::string bifParamName = pitchforkParams->get(
                          "Bifurcation Parameter",
                          "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamID[0] = p.getIndex(bifParamName);

  if (!pitchforkParams->isParameter("Antisymmetric Vector")) {
    globalData->locaErrorCheck->throwError(func,
               "\"Antisymmetric Vector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> asymVecPtr =
    (*pitchforkParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<NOX::Abstract::Vector> >("Antisymmetric Vector");

  if (!pitchforkParams->isParameter("Length Normalization Vector")) {
    globalData->locaErrorCheck->throwError(func,
               "\"Length Normalization Vector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> lenVecPtr =
    (*pitchforkParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<NOX::Abstract::Vector> >("Length Normalization Vector");

  if (!pitchforkParams->isParameter("Initial Null Vector")) {
    globalData->locaErrorCheck->throwError(func,
                 "\"Initial Null Vector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> nullVecPtr =
    (*pitchforkParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<NOX::Abstract::Vector> >("Initial Null Vector");

  bool perturbSoln = pitchforkParams->get(
                           "Perturb Initial Solution",
                           false);
  double perturbSize = pitchforkParams->get(
                         "Relative Perturbation Size",
                         1.0e-3);
  asymMultiVec =
    asymVecPtr->createMultiVector(1, NOX::DeepCopy);
  lengthMultiVec =
    lenVecPtr->createMultiVector(1, NOX::DeepCopy);
  *(xMultiVec.getColumn(0)->getNullVec()) = *nullVecPtr;

  // Instantiate solver strategy
  solverStrategy =
    globalData->locaFactory->createMooreSpencePitchforkSolverStrategy(
                                  parsedParams,
                                  pitchforkParams);

  // Set up multi-vector views
  setupViews();

  init(perturbSoln, perturbSize);
}

LOCA::Pitchfork::MooreSpence::ExtendedGroup::ExtendedGroup(
        const LOCA::Pitchfork::MooreSpence::ExtendedGroup& source,
        NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    pitchforkParams(source.pitchforkParams),
    grpPtr(Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MooreSpence::AbstractGroup>(source.grpPtr->clone(type))),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    asymMultiVec(source.asymMultiVec->clone(type)),
    lengthMultiVec(source.lengthMultiVec->clone(type)),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    asymVec(),
    lengthVec(),
    solverStrategy(source.solverStrategy),
    index_f(1),
    index_dfdp(1),
    bifParamID(source.bifParamID),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton)
{

  // Instantiate solver strategy
  solverStrategy =
    globalData->locaFactory->createMooreSpencePitchforkSolverStrategy(
                                     parsedParams,
                                 pitchforkParams);

  // Set up multi-vector views
  setupViews();

  if (type == NOX::ShapeCopy) {
    isValidF = false;
    isValidJacobian = false;
    isValidNewton = false;
  }
}

LOCA::Pitchfork::MooreSpence::ExtendedGroup::~ExtendedGroup()
{
}

NOX::Abstract::Group&
LOCA::Pitchfork::MooreSpence::ExtendedGroup::operator=(
                       const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Pitchfork::MooreSpence::ExtendedGroup::clone(NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedGroup(*this,
                                 type));
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setX(
                          const NOX::Abstract::Vector& y)
{
  const LOCA::Pitchfork::MooreSpence::ExtendedVector& yy =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedVector&>(y);
  grpPtr->setX( *yy.getXVec() );
  *xVec = y;
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeX(
                          const NOX::Abstract::Group& g,
                          const NOX::Abstract::Vector& d,
                          double step)
{
  const LOCA::Pitchfork::MooreSpence::ExtendedGroup& gg =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedGroup&>(g);
  const LOCA::Pitchfork::MooreSpence::ExtendedVector& dd =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedVector&>(d);

  grpPtr->computeX(*(gg.grpPtr), *dd.getXVec(), step);
  xVec->update(1.0, gg.getX(), step, dd, 0.0);
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeF()
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  double sigma = xVec->getSlack();

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }
  fVec->getXVec()->update(1.0, grpPtr->getF(), sigma, *asymVec, 0.0);

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Compute J*n
  status = grpPtr->applyJacobian(*(xVec->getNullVec()),
                 *(fVec->getNullVec()));
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  // Compute <x,psi>
  fVec->getSlack() = grpPtr->innerProduct(*(xVec->getXVec()), *asymVec);

  // Compute phi^T*n
  fVec->getBifParam() = lTransNorm(*(xVec->getNullVec())) - 1.0;

  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeJacobian()
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  // Note:  the first column of fMultiVec stores f + sigma*psi, not f,
  // so we always need to recompute f.  This changes the first column
  // of fMultiVec back to f
  status = grpPtr->computeDfDpMulti(bifParamID,
                    *fMultiVec.getXMultiVec(),
                    false);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  // Change the first column back to f + sigma*psi
  double sigma = xVec->getSlack();
  fVec->getXVec()->update(sigma, *asymVec, 1.0);

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
          asymMultiVec,
          xVec->getNullVec(),
          fVec->getNullVec(),
          fMultiVec.getColumn(1)->getXVec(),
          fMultiVec.getColumn(1)->getNullVec());

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeGradient()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeNewton(
                         Teuchos::ParameterList& params)
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeNewton()";
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
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobian(
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
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianTranspose(
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
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianInverse(
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
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianMultiVector(
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
                 "Called with invalid Jacobian!");
  }

  // Cast vectors to pitchfork vectors
  const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& pf_input =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& pf_result =
    dynamic_cast<LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&>(result);

  // Get constant references to input vector components
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x =
    pf_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_null =
    pf_input.getNullMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_slack = pf_input.getSlacks();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = pf_input.getBifParams();

  // Get non-constant references to result vector components
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x =
    pf_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_null =
    pf_result.getNullMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_slack =
    pf_result.getSlacks();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param =
    pf_result.getBifParams();

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

  // compute J*x + psi*sigma
  result_x->update(Teuchos::NO_TRANS, 1.0, *asymMultiVec, *input_slack);


  // compute J*x + psi*sigma + (dR/dp)*p
  result_x->update(Teuchos::NO_TRANS, 1.0, *(dfdpMultiVec->getXMultiVec()),
           *input_param);

  // compute J*y
  status = grpPtr->applyJacobianMultiVector(*input_null, *result_null);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // compute J*y + (dJy/dp)*p
  result_null->update(Teuchos::NO_TRANS, 1.0,
              *(dfdpMultiVec->getNullMultiVec()),
              *input_param);

  // compute (dJy/dx)*x
  status = grpPtr->computeDJnDxaMulti(*(xVec->getNullVec()),
                      *(fVec->getNullVec()),
                      *input_x, *tmp);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null->update(1.0, *tmp, 1.0);

  // compute <x,psi>
  grpPtr->innerProduct(*asymMultiVec, *input_x, *result_slack);

  // compute l^T*y
  lTransNorm(*input_null, *result_param);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector(
                     const NOX::Abstract::MultiVector& /* input */,
                     NOX::Abstract::MultiVector& /* result */) const
{
  globalData->locaErrorCheck->throwError(
          "LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector()",
          "Method not implemented!");

  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector(
                     Teuchos::ParameterList& params,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  // Cast vectors to TP vectors
  const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& pf_input =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& pf_result =
    dynamic_cast<LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&>(result);

  NOX::Abstract::Group::ReturnType res = solverStrategy->solve(params, pf_input, pf_result);

  return res;
}

bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isF() const
{
  return isValidF;
}

bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isJacobian() const
{
  return isValidJacobian;
}

bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isGradient() const
{
  return false;
}

bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isNewton() const
{
  return isValidNewton;
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getX() const
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getF() const
{
  return *fVec;
}

double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNormF() const
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradient() const
{
  globalData->locaErrorCheck->throwError(
          "LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradient()",
          " - not implemented");
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNewton() const
{
  return *newtonVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getXPtr() const
{
  return xVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getFPtr() const
{
  return fVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradientPtr() const
{
  globalData->locaErrorCheck->throwError(
          "LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradientPtr()",
          " - not implemented");
  return getNewtonPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNewtonPtr() const
{
  return newtonVec;
}

double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual() const
{
  std::string callingFunction =
    "LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::Pitchfork::MooreSpence::ExtendedVector residual = *fVec;

  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getUnderlyingGroup()
{
  return grpPtr;
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::copy(
                        const NOX::Abstract::Group& src)
{
  const LOCA::Pitchfork::MooreSpence::ExtendedGroup& source =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedGroup&>(src);

  // Protect against A = A
  if (this != &source) {

    // Copy values
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    pitchforkParams = source.pitchforkParams;
    grpPtr->copy(*(source.grpPtr));
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    *asymMultiVec = *source.asymMultiVec;
    *lengthMultiVec = *source.lengthMultiVec;
    index_f = source.index_f;
    index_dfdp = source.index_dfdp;
    bifParamID = source.bifParamID;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;

    // set up views again just to be safe
    setupViews();

    // Instantiate solver strategy
    solverStrategy =
      globalData->locaFactory->createMooreSpencePitchforkSolverStrategy(
                                 parsedParams,
                                             pitchforkParams);
  }
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParamsMulti(
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
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParams(
                          const LOCA::ParameterVector& p)
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
  setBifParam(p[bifParamID[0]]);
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParam(int paramID,
                              double val)
{
  if (paramID == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParam(std::string paramID,
                              double val)
{
  const LOCA::ParameterVector& pVec = grpPtr->getParams();
  if (pVec.getIndex(paramID) == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

const LOCA::ParameterVector&
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParams() const
{
  return grpPtr->getParams();
}

double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParam(std::string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeDfDpMulti(
                        const std::vector<int>& paramIDs,
                        NOX::Abstract::MultiVector& dfdp,
                        bool isValid_F)
{
   std::string callingFunction =
    "LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeDfDpMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast dfdp to pitchfork vector
  LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& pf_dfdp =
    dynamic_cast<LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&>(dfdp);

  // Compute df/dp
  // Note:  the first column of fMultiVec stores f + sigma*psi, not f,
  // so we always need to recompute f.  This changes the first column
  // of fMultiVec back to f
  status = grpPtr->computeDfDpMulti(paramIDs, *pf_dfdp.getXMultiVec(),
                    false);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Change the first column back to f + sigma*psi
  double sigma = xVec->getSlack();
  pf_dfdp.getColumn(0)->getXVec()->update(sigma, *asymVec, 1.0);

  // Compute d(Jn)/dp
  status = grpPtr->computeDJnDpMulti(paramIDs, *(xVec->getNullVec()),
                     *pf_dfdp.getNullMultiVec(), isValid_F);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Set parameter components
  if (!isValid_F) {
    pf_dfdp.getScalar(0,0) = grpPtr->innerProduct(*(xVec->getXVec()),
                          *asymVec);
    pf_dfdp.getScalar(1,0) = lTransNorm(*(xVec->getNullVec()));
  }
  for (int i=0; i<dfdp.numVectors()-1; i++) {
    pf_dfdp.getScalar(0,i+1) = 0.0;
    pf_dfdp.getScalar(1,i+1) = 0.0;
  }

  return finalStatus;
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::preProcessContinuationStep(
             LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->preProcessContinuationStep(stepStatus);
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::postProcessContinuationStep(
             LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->postProcessContinuationStep(stepStatus);
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::projectToDraw(
                           const NOX::Abstract::Vector& x,
                           double *px) const
{
  const LOCA::Pitchfork::MooreSpence::ExtendedVector& mx =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedVector&>(x);

  grpPtr->projectToDraw(*(mx.getXVec()), px);
  px[grpPtr->projectToDrawDimension()] = mx.getBifParam();
}

int
LOCA::Pitchfork::MooreSpence::ExtendedGroup::projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension() + 1;
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution(
                          const double conParam) const
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Pitchfork located at: " <<
      globalData->locaUtils->sciformat(conParam) << "   " <<
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;

    globalData->locaUtils->out() <<
      "\tSlack variable sigma = " <<
      globalData->locaUtils->sciformat(xVec->getSlack()) << std::endl;

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
LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution(
                         const NOX::Abstract::Vector& x_,
                         const double conParam) const
{
  const LOCA::Pitchfork::MooreSpence::ExtendedVector& pf_x =
    dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedVector&>(x_);

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Pitchfork located at: " <<
      globalData->locaUtils->sciformat(conParam) << "   " <<
      globalData->locaUtils->sciformat(pf_x.getBifParam()) << std::endl;

    globalData->locaUtils->out() <<
      "\tSlack variable sigma = " <<
      globalData->locaUtils->sciformat(pf_x.getSlack()) << std::endl;

    globalData->locaUtils->out() <<
      "\tPrinting Solution Vector for conParam = " <<
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(*pf_x.getXVec(), conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tPrinting Null Vector for bif param = " <<
      globalData->locaUtils->sciformat(pf_x.getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*pf_x.getNullVec(), pf_x.getBifParam());
}

double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getBifParam() const
{
  return grpPtr->getParam(bifParamID[0]);
}

double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::lTransNorm(
                    const NOX::Abstract::Vector& n) const
{
  return lengthVec->innerProduct(n) / lengthVec->length();
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::lTransNorm(
            const NOX::Abstract::MultiVector& n,
            NOX::Abstract::MultiVector::DenseMatrix& result) const
{
  n.multiply(1.0 / lengthVec->length(), *lengthMultiVec, result);
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setBifParam(double param)
{
  grpPtr->setParam(bifParamID[0], param);
  xVec->getBifParam() = param;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setupViews()
{
  index_f[0] = 0;
  index_dfdp[0] = 1;

  xVec = xMultiVec.getColumn(0);
  fVec = fMultiVec.getColumn(0);
  newtonVec = newtonMultiVec.getColumn(0);
  asymVec = Teuchos::rcp(&(*asymMultiVec)[0], false);
  lengthVec = Teuchos::rcp(&(*lengthMultiVec)[0],false);

  ffMultiVec = Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MooreSpence::ExtendedMultiVector>(fMultiVec.subView(index_f),true);

  dfdpMultiVec = Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MooreSpence::ExtendedMultiVector>(fMultiVec.subView(index_dfdp),true);

}

void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::init(bool perturbSoln,
                          double perturbSize)
{
  xVec->getBifParam() = getBifParam();

  // Rescale length vector so that the normalization condition is met
  double lVecDotNullVec = lTransNorm(*(xVec->getNullVec()));

  if (fabs(lVecDotNullVec) < 1.0e-8) {
    globalData->locaErrorCheck->throwError(
           "LOCA::Pitchfork::MooreSpence::ExtendedGroup::init()",
           "null vector cannot be orthogonal to length-scaling vector: ");
  }
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tIn LOCA::Pitchfork::MooreSpence::ExtendedGroup::init(), " <<
      "scaling null vector by:" <<
      globalData->locaUtils->sciformat(1.0 / lVecDotNullVec) << std::endl;
  }
  xVec->getNullVec()->scale(1.0/lVecDotNullVec);

  // Rescale asymmetric vector to have unit length
  double psi_norm = sqrt( grpPtr->innerProduct(*asymVec, *asymVec) );
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() <<
      "\tIn LOCA::Pitchfork::MooreSpence::ExtendedGroup::init(), " <<
      "scaling asymmetric vector by:" <<
      globalData->locaUtils->sciformat(1.0 / psi_norm) << std::endl;
  }
  asymVec->scale(1.0 / psi_norm);

  if (perturbSoln) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
     globalData->locaUtils->out() <<
       "\tIn LOCA::Pitchfork::MooreSpence::ExtendedGroup::init(), " <<
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

