// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_TurningPoint_MinimallyAugmented_Constraint.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_BorderedSolver_JacobianOperator.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"
#include "LOCA_TimeDependent_AbstractGroup.H"

LOCA::TurningPoint::MinimallyAugmented::Constraint::
Constraint(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
    const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g,
    int bif_param) :
  globalData(global_data),
  parsedParams(topParams),
  turningPointParams(tpParams),
  grpPtr(g),
  a_vector(grpPtr->getX().createMultiVector(1, NOX::ShapeCopy)),
  b_vector(a_vector->clone(NOX::ShapeCopy)),
  w_vector(a_vector->clone(NOX::ShapeCopy)),
  v_vector(a_vector->clone(NOX::ShapeCopy)),
  Jv_vector(a_vector->clone(NOX::ShapeCopy)),
  Jtw_vector(a_vector->clone(NOX::ShapeCopy)),
  sigma_x(a_vector->clone(NOX::ShapeCopy)),
  constraints(1, 1),
  borderedSolver(),
  dn(static_cast<double>(a_vector->length())),
  sigma_scale(1.0),
  isSymmetric(false),
  isValidConstraints(false),
  isValidDX(false),
  bifParamID(1),
  updateVectorsEveryContinuationStep(true),
  updateVectorsEveryIteration(false),
  nullVecScaling(NVS_OrderN),
  multiplyMass(false),
  tdGrp(Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grpPtr)),
  tmp_mass(grpPtr->getX().clone(NOX::ShapeCopy))
{
  // Instantiate bordered solvers
  borderedSolver =
    globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
                              turningPointParams);

  // Get symmetric flag
  isSymmetric = turningPointParams->get("Symmetric Jacobian", false);

  bifParamID[0] = bif_param;

  // Options
  updateVectorsEveryContinuationStep =
    turningPointParams->get("Update Null Vectors Every Continuation Step",
                true);
  updateVectorsEveryIteration =
    turningPointParams->get("Update Null Vectors Every Nonlinear Iteration",
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
       "LOCA::TurningPoint::MinimallyAugmented::Constraint::Constraint()",
       std::string("Unknown null vector scaling method:  ") + nullVecScalingMethod);
  multiplyMass =
    turningPointParams->get("Multiply Null Vectors by Mass Matrix", false);
  if (multiplyMass && tdGrp == Teuchos::null) {
    globalData->locaErrorCheck->throwError(
       "LOCA::TurningPoint::MinimallyAugmented::Constraint::Constraint()",
       "Group must be derived from LOCA::TimeDependent::AbstractGroup to multiply null vectors by mass matrix");
  }

  // Compute/get initial "a" & "b" vectors
  getInitialVectors((*a_vector)[0], (*b_vector)[0]);
}

LOCA::TurningPoint::MinimallyAugmented::Constraint::
Constraint(const LOCA::TurningPoint::MinimallyAugmented::Constraint& source,
       NOX::CopyType type) :
  globalData(source.globalData),
  parsedParams(source.parsedParams),
  turningPointParams(source.turningPointParams),
  grpPtr(Teuchos::null),
  a_vector(source.a_vector->clone(type)),
  b_vector(source.b_vector->clone(type)),
  w_vector(source.w_vector->clone(type)),
  v_vector(source.v_vector->clone(type)),
  Jv_vector(source.Jv_vector->clone(type)),
  Jtw_vector(source.Jtw_vector->clone(type)),
  sigma_x(source.sigma_x->clone(type)),
  constraints(source.constraints),
  borderedSolver(),
  dn(source.dn),
  sigma_scale(source.sigma_scale),
  isSymmetric(source.isSymmetric),
  isValidConstraints(false),
  isValidDX(false),
  bifParamID(source.bifParamID),
  updateVectorsEveryContinuationStep(source.updateVectorsEveryContinuationStep),
  updateVectorsEveryIteration(source.updateVectorsEveryIteration),
  nullVecScaling(source.nullVecScaling),
  multiplyMass(source.multiplyMass),
  tdGrp(source.tdGrp),
  tmp_mass(source.tmp_mass->clone(type))
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
  if (source.isValidDX && type == NOX::DeepCopy)
    isValidDX = true;

  // Instantiate bordered solvers
  borderedSolver =
    globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
                              turningPointParams);

  // We don't explicitly copy the group because the constrained group
  // will do that
}

LOCA::TurningPoint::MinimallyAugmented::Constraint::
~Constraint()
{
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setGroup(const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g)
{
  grpPtr = g;
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getLeftNullVec() const
{
  return Teuchos::rcp(&(*w_vector)[0], false);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getRightNullVec() const
{
  return Teuchos::rcp(&(*v_vector)[0], false);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getAVec() const
{
  return Teuchos::rcp(&(*a_vector)[0], false);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getBVec() const
{
  return Teuchos::rcp(&(*b_vector)[0], false);
}

double
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getSigma() const
{
  return constraints(0,0);
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::TurningPoint::MinimallyAugmented::Constraint& source =
  dynamic_cast<const LOCA::TurningPoint::MinimallyAugmented::Constraint&>(src);

  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    turningPointParams = source.turningPointParams;
    *a_vector = *(source.a_vector);
    *b_vector = *(source.b_vector);
    *w_vector = *(source.w_vector);
    *v_vector = *(source.v_vector);
    *Jv_vector = *(source.Jv_vector);
    *Jtw_vector = *(source.Jtw_vector);
    *sigma_x = *(source.sigma_x);
    constraints.assign(source.constraints);
    dn = source.dn;
    sigma_scale = source.sigma_scale;
    isSymmetric = source.isSymmetric;
    isValidConstraints = source.isValidConstraints;
    isValidDX = source.isValidDX;
    bifParamID = source.bifParamID;
    updateVectorsEveryContinuationStep =
      source.updateVectorsEveryContinuationStep;
    updateVectorsEveryIteration =
      source.updateVectorsEveryIteration;
    nullVecScaling = source.nullVecScaling;
    multiplyMass = source.multiplyMass;

    // Instantiate bordered solvers
    borderedSolver =
      globalData->locaFactory->createBorderedSolverStrategy(
                             parsedParams,
                             turningPointParams);

    // We don't explicitly copy the group because the constrained group
    // will do that
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Constraint(*this, type));
}

int
LOCA::TurningPoint::MinimallyAugmented::Constraint::
numConstraints() const
{
  return 1;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setX(const NOX::Abstract::Vector& y)
{
  grpPtr->setX(y);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setParams(const std::vector<int>& paramIDs,
      const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  isValidConstraints = false;
  isValidDX = false;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::
computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute J
  status = grpPtr->computeJacobian();
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Set up bordered systems
  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> op =
    Teuchos::rcp(new  LOCA::BorderedSolver::JacobianOperator(grpPtr));
  borderedSolver->setMatrixBlocksMultiVecConstraint(op,
                            a_vector,
                            b_vector,
                            Teuchos::null);

  // Create RHS
  NOX::Abstract::MultiVector::DenseMatrix one(1,1);
  if (nullVecScaling == NVS_OrderN)
    one(0,0) = dn;
  else
    one(0,0) = 1.0;

  // Get linear solver parameters
  Teuchos::RCP<Teuchos::ParameterList> linear_solver_params =
    parsedParams->getSublist("Linear Solver");

  // Compute sigma_1 and right null vector v
  NOX::Abstract::MultiVector::DenseMatrix s1(1,1);
  status = borderedSolver->initForSolve();
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);
  status = borderedSolver->applyInverse(*linear_solver_params,
                    NULL,
                    &one,
                    *v_vector,
                    s1);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Compute sigma_2 and left null vector w
  NOX::Abstract::MultiVector::DenseMatrix s2(1,1);
  if (!isSymmetric) {
    status = borderedSolver->initForTransposeSolve();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
    status = borderedSolver->applyInverseTranspose(*linear_solver_params,
                           NULL,
                           &one,
                           *w_vector,
                           s2);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

  }
  else {
    *w_vector = *v_vector;
    s2.assign(s1);
  }

  // Compute sigma = -w^T*J*v
  status = grpPtr->applyJacobianMultiVector(*v_vector, *Jv_vector);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);
  if (!isSymmetric) {
    status = grpPtr->applyJacobianTransposeMultiVector(*w_vector, *Jtw_vector);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }
  else
    *Jtw_vector = *Jv_vector;
  Jv_vector->multiply(-1.0, *w_vector, constraints);

  // Scale sigma
  double w_norm = (*w_vector)[0].norm();
  double v_norm = (*v_vector)[0].norm();
  double Jv_norm = (*Jv_vector)[0].norm();
  double Jtw_norm = (*Jtw_vector)[0].norm();
  if (nullVecScaling == NVS_OrderN)
    sigma_scale = dn;
  else
    sigma_scale = 1.0;
  constraints.scale(1.0/sigma_scale);

  if (globalData->locaUtils->isPrintType(NOX::Utils::OuterIteration)) {
    globalData->locaUtils->out() << "\n\t||Right null vector v|| = "
                 << globalData->locaUtils->sciformat(v_norm);
    globalData->locaUtils->out() << "\n\t||Left null vector w|| = "
                 << globalData->locaUtils->sciformat(w_norm);
    globalData->locaUtils->out() << "\n\t||Jv|| = "
                 << globalData->locaUtils->sciformat(Jv_norm);
    globalData->locaUtils->out() << "\n\t||J^T*w|| = "
                 << globalData->locaUtils->sciformat(Jtw_norm);
    globalData->locaUtils->out() <<
      "\n\tRight estimate for singularity of Jacobian (sigma1) = " <<
      globalData->locaUtils->sciformat(s1(0,0));
    globalData->locaUtils->out() <<
      "\n\tLeft estimate for singularity of Jacobian (sigma2) = " <<
      globalData->locaUtils->sciformat(s2(0,0));
    globalData->locaUtils->out() <<
      "\n\tFinal Estimate for singularity of Jacobian (sigma) = " <<
      globalData->locaUtils->sciformat(constraints(0,0)) << std::endl;
  }

  isValidConstraints = true;

  // Update a and b if requested
  if (updateVectorsEveryIteration) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::OuterIteration)) {
      globalData->locaUtils->out() <<
    "\n\tUpdating null vectors for the next nonlinear iteration" <<
    std::endl;
    }
    *a_vector = *w_vector;
    *b_vector = *v_vector;

    scaleNullVectors((*a_vector)[0],(*b_vector)[0]);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::
computeDX()
{
  if (isValidDX)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDX()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute sigma, w and v if necessary
  if (!isValidConstraints) {
    status = computeConstraints();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Compute -(w^T*J*v)_x
  status = grpPtr->computeDwtJnDx((*w_vector)[0], (*v_vector)[0],
                  (*sigma_x)[0]);
  finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  sigma_x->scale(-1.0/sigma_scale);

  isValidDX = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::
computeDP(const std::vector<int>& paramIDs,
      NOX::Abstract::MultiVector::DenseMatrix& dgdp,
      bool /* isValidG */)
{
  std::string callingFunction =
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDP()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute sigma, w and v if necessary
  if (!isValidConstraints) {
    status = computeConstraints();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Compute -(w^T*J*v)_p
  status = grpPtr->computeDwtJnDp(paramIDs, (*w_vector)[0], (*v_vector)[0],
                  dgdp, false);
  finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  dgdp.scale(-1.0/sigma_scale);

  // Set the first column of dgdp
  dgdp(0,0) = constraints(0,0);

  return finalStatus;
}

bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::
isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::
isDX() const
{
  return isValidDX;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getDX() const
{
  return sigma_x.get();
}

bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::
isDXZero() const
{
  return false;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful &&
      updateVectorsEveryContinuationStep) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
      globalData->locaUtils->out() <<
      "\n\tUpdating null vectors for the next continuation step" << std::endl;
    }
    *a_vector = *w_vector;
    *b_vector = *v_vector;

    scaleNullVectors((*a_vector)[0],(*b_vector)[0]);
  }
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
scaleNullVectors(NOX::Abstract::Vector& a, NOX::Abstract::Vector& b)
{
  std::string callingFunction =
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::scaleNullVectors()";
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
    *tmp_mass = b;
    status = tdGrp->applyShiftedMatrix(*tmp_mass, b);
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
    b.scale(1.0 / b.norm());
  }
  else if (nullVecScaling == NVS_OrderN) {
    a.scale(std::sqrt(dn) / a.norm());
    b.scale(std::sqrt(dn) / b.norm());
  }
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getInitialVectors(NOX::Abstract::Vector& aVec,
          NOX::Abstract::Vector& bVec)
{
  std::string callingFunction =
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::getIntitialVectors()";

  // Get method
  std::string method =
    turningPointParams->get("Initial Null Vector Computation",
                "User Provided");
  if (method == "Solve df/dp") {
    NOX::Abstract::Group::ReturnType status;
    NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
    std::vector<int> paramID(1);
    paramID[0] = bifParamID[0];
    Teuchos::RCP<NOX::Abstract::MultiVector> fdfdp =
      grpPtr->getX().createMultiVector(2);
    aVec.init(0.0);
    bVec.init(0.0);

    // Compute df/dp
    status = grpPtr->computeDfDpMulti(paramID, *fdfdp, false);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Compute J
    status = grpPtr->computeJacobian();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Compute b = J^-1*dfdp
    Teuchos::RCP<Teuchos::ParameterList> lsParams =
      parsedParams->getSublist("Linear Solver");
    status = grpPtr->applyJacobianInverse(*lsParams, (*fdfdp)[1], bVec);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Compute a = J^-T*dfdp if necessary
    if (!isSymmetric) {
      // Cast group to one that can solve J^T
      Teuchos::RCP<LOCA::Abstract::TransposeSolveGroup> ts_grp =
    Teuchos::rcp_dynamic_cast<LOCA::Abstract::TransposeSolveGroup>(grpPtr);
      if (ts_grp == Teuchos::null)
    globalData->locaErrorCheck->throwError(
       callingFunction,
       std::string("Group must implement LOCA::Abstract::TransposeSolveGroup") +
       std::string(" to compute initial left null vector"));

      Teuchos::RCP<Teuchos::ParameterList> lsParams2 =
    parsedParams->getSublist("Linear Solver");
      status =
    ts_grp->applyJacobianTransposeInverse(*lsParams2, (*fdfdp)[1], aVec);
      finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(
                                 status,
                                 finalStatus,
                                 callingFunction);
    }
    else
      aVec = bVec;
  }
  else if (method == "Constant") {
    aVec.init(1.0);
    bVec.init(1.0);
  }

  else {

    // Get initial "a" vector
    if (!turningPointParams->isParameter("Initial A Vector")) {
      globalData->locaErrorCheck->throwError(callingFunction,
                     "\"Initial A Vector\" is not set!");
    }
    Teuchos::RCP<NOX::Abstract::Vector> aVecPtr =
      (*turningPointParams).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<NOX::Abstract::Vector> >("Initial A Vector");
    aVec = *aVecPtr;

    // Get initial "b" vector
    if (!isSymmetric) {
      if (!turningPointParams->isParameter("Initial B Vector")) {
    globalData->locaErrorCheck->throwError(callingFunction,
                       "\"Initial B Vector\" is not set!");
      }
      Teuchos::RCP<NOX::Abstract::Vector> bVecPtr =
    (*turningPointParams).INVALID_TEMPLATE_QUALIFIER
        get< Teuchos::RCP<NOX::Abstract::Vector> >("Initial B Vector");
      bVec = *bVecPtr;
    }
    else
      bVec = aVec;
  }

  scaleNullVectors(aVec, bVec);
}
