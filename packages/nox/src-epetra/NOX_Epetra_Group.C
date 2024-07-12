// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Epetra_Group.H"    // class definition

#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_Epetra_Interface_Required.H"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "AztecOO_ConditionNumber.h"
#include "NOX_SolverStats.hpp"

using namespace NOX;
using namespace NOX::Epetra;
using namespace Teuchos;

Group::Group(Teuchos::ParameterList& printParams,
         const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
         const NOX::Epetra::Vector& x):
  utils(printParams),
  xVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(DeepCopy))),
  xVector(*xVectorPtr),
  RHSVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(ShapeCopy))),
  RHSVector(*RHSVectorPtr),
  gradVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(ShapeCopy))),
  gradVector(*gradVectorPtr),
  NewtonVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(ShapeCopy))),
  NewtonVector(*NewtonVectorPtr),
  normNewtonSolveResidual(0),
  conditionNumber(0.0),
  sharedLinearSystem(*sharedLinearSystemPtr),
  linearResidCompDisabled(false),
  userInterfacePtr(i),
  linearSolveConverged(false),
  numIterations(-1),
  achievedTol(-1.0)
{
  // Set all isValid flags to false
  resetIsValid();
}

Group::Group(Teuchos::ParameterList& printParams,
         const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
         const NOX::Epetra::Vector& x,
         const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys):
  utils(printParams),
  xVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(DeepCopy))),
  xVector(*xVectorPtr),
  RHSVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(ShapeCopy))),
  RHSVector(*RHSVectorPtr),
  gradVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(ShapeCopy))),
  gradVector(*gradVectorPtr),
  NewtonVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(x.clone(ShapeCopy))),
  NewtonVector(*NewtonVectorPtr),
  normNewtonSolveResidual(0),
  conditionNumber(0.0),
  sharedLinearSystemPtr(Teuchos::rcp(new SharedObject<NOX::Epetra::LinearSystem, NOX::Epetra::Group>(linSys))),
  sharedLinearSystem(*sharedLinearSystemPtr),
  linearResidCompDisabled(false),
  userInterfacePtr(i),
  linearSolveConverged(false),
  numIterations(-1),
  achievedTol(-1.0)
{
  // Set all isValid flags to false
  resetIsValid();
}

Group::Group(const Group& source, CopyType type) :
  utils(source.utils),
  xVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(source.xVector.clone(type))),
  xVector(*xVectorPtr),
  RHSVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(source.RHSVector.clone(type))),
  RHSVector(*RHSVectorPtr),
  gradVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(source.gradVector.clone(type))),
  gradVector(*gradVectorPtr),
  NewtonVectorPtr(rcp_dynamic_cast<NOX::Epetra::Vector>(source.NewtonVector.clone(type))),
  NewtonVector(*NewtonVectorPtr),
  normNewtonSolveResidual(0.0),
  conditionNumber(0.0),
  sharedLinearSystemPtr(source.sharedLinearSystemPtr),
  sharedLinearSystem(*sharedLinearSystemPtr),
  linearResidCompDisabled(source.linearResidCompDisabled),
  userInterfacePtr(source.userInterfacePtr),
  linearSolveConverged(source.linearSolveConverged),
  numIterations(source.numIterations),
  achievedTol(source.achievedTol)
{

  switch (type) {

  case DeepCopy:

    isValidRHS = source.isValidRHS;
    isValidJacobian = source.isValidJacobian;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
    isValidConditionNumber = source.isValidConditionNumber;
    normNewtonSolveResidual = source.normNewtonSolveResidual;
    conditionNumber = source.conditionNumber;
    isValidPreconditioner = source.isValidPreconditioner;
    isValidSolverJacOp = source.isValidSolverJacOp;

    // New copy takes ownership of the shared Jacobian for DeepCopy
    if (isValidJacobian)
      sharedLinearSystem.getObject(this);

    break;

  case ShapeCopy:
    resetIsValid();
    break;

  default:
    std::cerr << "ERROR: Invalid ConstructorType for group copy constructor." << std::endl;
    throw std::runtime_error("NOX Error");
  }

}

Group::~Group()
{
}

void Group::resetIsValid() //private
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
  isValidNormNewtonSolveResidual = false;
  isValidPreconditioner = false;
  isValidSolverJacOp = false;
  isValidConditionNumber = false;
  return;
}

Teuchos::RCP<NOX::Abstract::Group> Group::clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp =
    Teuchos::rcp(new NOX::Epetra::Group(*this, type));
  return newgrp;
}

Abstract::Group& Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group& Group::operator=(const Group& source)
{
  // Copy the xVector
  xVector = source.xVector;

  // Copy reference to sharedJacobian
  //sharedLinearSystemPtr = source.sharedLinearSystemPtr;

  // Update the isValidVectors
  isValidRHS = source.isValidRHS;
  isValidGrad = source.isValidGrad;
  isValidNewton = source.isValidNewton;
  isValidJacobian = source.isValidJacobian;
  isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
  isValidPreconditioner = source.isValidPreconditioner;
  isValidSolverJacOp = source.isValidSolverJacOp;
  isValidConditionNumber = source.isValidConditionNumber;

  // Only copy vectors that are valid
  if (isValidRHS) {
    RHSVector = source.RHSVector;
  }

  if (isValidGrad)
    gradVector = source.gradVector;

  if (isValidNewton)
    NewtonVector = source.NewtonVector;

  if (isValidNormNewtonSolveResidual)
    normNewtonSolveResidual = source.normNewtonSolveResidual;

  // If valid, this takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedLinearSystem.getObject(this);


  if (isValidConditionNumber)
    conditionNumber = source.conditionNumber;

  linearResidCompDisabled = source.linearResidCompDisabled;

  linearSolveConverged = source.linearSolveConverged;
  numIterations = source.numIterations;
  achievedTol = source.achievedTol;

  return *this;
}

void Group::setX(const Abstract::Vector& y)
{
  setX(dynamic_cast<const NOX::Epetra::Vector&> (y));
  return;
}

void Group::setX(const NOX::Epetra::Vector& y)
{
  if (isPreconditioner()) {
    sharedLinearSystem.getObject(this)->destroyPreconditioner();
  }
  resetIsValid();
  xVector = y;
  return;
}

void Group::computeX(const NOX::Abstract::Group& grp,
             const NOX::Abstract::Vector& d,
             double step)
{
  // Cast to appropriate type, then call the "native" computeX
  const NOX::Epetra::Group& epetragrp = dynamic_cast<const Group&> (grp);
  const NOX::Epetra::Vector& epetrad =
    dynamic_cast<const NOX::Epetra::Vector&> (d);
  computeX(epetragrp, epetrad, step);
  return;
}

void Group::computeX(const Group& grp,
             const NOX::Epetra::Vector& d,
             double step)
{
  if (isPreconditioner())
    sharedLinearSystem.getObject(this)->destroyPreconditioner();
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return;
}

Abstract::Group::ReturnType Group::computeF()
{
  if (isF())
    return Abstract::Group::Ok;

  bool status = false;

  status = userInterfacePtr->computeF(xVector.getEpetraVector(), RHSVector.getEpetraVector(), NOX::Epetra::Interface::Required::Residual);

  if (status == false) {
    std::cout << "ERROR: Epetra::Group::computeF() - fill failed!!!"
     << std::endl;
    throw std::runtime_error("NOX Error: Fill Failed");
  }

  isValidRHS = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeJacobian()
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return Abstract::Group::Ok;

  // Fill the Jacobian
  bool status = false;

  status = sharedLinearSystem.getObject(this)->
    computeJacobian(xVector);

  if (status == false) {
    std::cout << "ERROR: NOX::Epetra::Group::computeJacobian() - fill failed!!!"
     << std::endl;
    throw std::runtime_error("NOX Error: Fill Failed");
  }

  // Update status of Jacobian wrt solution vector
  isValidJacobian = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeGradient()
{
  if (isGradient())
    return Abstract::Group::Ok;

  if (!isF()) {
    std::cerr << "ERROR: NOX::Epetra::Group::computeGradient() - RHS is out of date wrt X!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::Epetra::Group::computeGradient() - Jacobian is out of date wrt X!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Compute grad = Jacobian^T * RHS.
  sharedLinearSystem.getObject(this)->applyJacobianTranspose(RHSVector,
                                 gradVector);

  // Update state
  isValidGrad = true;

  // Return result
  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeNewton(Teuchos::ParameterList& p)
{
  if (isNewton())
    return Abstract::Group::Ok;

  if (!isF()) {
    std::cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid RHS" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid Jacobian" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  Abstract::Group::ReturnType status;

  // Zero out the Newton Vector
  NewtonVector.init(0.0);

  // Create Epetra problem for the linear solve
  status = applyJacobianInverse(p, RHSVector, NewtonVector);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state EVEN IF LINEAR SOLVE FAILED
  // We still may want to use the vector even it it just missed it's
  isValidNewton = true;

  // Compute the 2-norm of the linear solve residual ||Js+f||
  // Can be disabled, but then disallows inexact Newton methods
  if (!linearResidCompDisabled)
    computeNormNewtonSolveResidual();

  // Return solution
  return status;
}

Abstract::Group::ReturnType Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetrainput =
    dynamic_cast<const NOX::Epetra::Vector&> (input);
  NOX::Epetra::Vector& epetraresult =
    dynamic_cast<NOX::Epetra::Vector&> (result);
  return applyJacobian(epetrainput, epetraresult);
}

Abstract::Group::ReturnType Group::applyJacobian(const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian())
    return Abstract::Group::BadDependency;

  // Apply the Jacobian
  bool status = sharedLinearSystem.getObject()->applyJacobian(input, result);

  return status == true ? Abstract::Group::Ok : Abstract::Group::Failed;
}

Abstract::Group::ReturnType Group::applyJacobianInverse (Teuchos::ParameterList &p, const Abstract::Vector &input, Abstract::Vector &result) const
{
  const NOX::Epetra::Vector& epetraInput = dynamic_cast<const NOX::Epetra::Vector&>(input);
  NOX::Epetra::Vector& epetraResult = dynamic_cast<NOX::Epetra::Vector&>(result);
  return applyJacobianInverse(p, epetraInput, epetraResult);
}

Abstract::Group::ReturnType Group::applyJacobianInverse (Teuchos::ParameterList &p, const NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const
{
  if (!isJacobian())
    return Abstract::Group::BadDependency;

  if (!isValidSolverJacOp) {
    sharedLinearSystem.getObject(this)->setJacobianOperatorForSolve(sharedLinearSystem.getObject(this)->getJacobianOperator());
    isValidSolverJacOp = true;
  }

  // Compute the preconditioner
  NOX::Epetra::LinearSystem::PreconditionerReusePolicyType precPolicy =
    sharedLinearSystem.getObject(this)->getPreconditionerPolicy();

  if (!isPreconditioner()) {
    if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REBUILD) {
      sharedLinearSystem.getObject(this)->destroyPreconditioner();
      sharedLinearSystem.getObject(this)->
    createPreconditioner(xVector, p, false);
      isValidPreconditioner = true;
    }
    else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_RECOMPUTE) {
      sharedLinearSystem.getObject(this)->recomputePreconditioner(xVector, p);
      isValidPreconditioner = true;
    }
    else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REUSE) {
      // Do Nothing!!!
    }
  }

  // Save linear solve stats
  linearSolveConverged = sharedLinearSystem.getObject(this)->applyJacobianInverse(p, input, result);
  numIterations = p.sublist("Output").get("Number of Linear Iterations",0);
  achievedTol = p.sublist("Output").get("Achieved Tolerance",0.0);

  return linearSolveConverged == true ? Abstract::Group::Ok : Abstract::Group::NotConverged;
}

Abstract::Group::ReturnType Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetrainput = dynamic_cast<const NOX::Epetra::Vector&> (input);
  NOX::Epetra::Vector& epetraresult = dynamic_cast<NOX::Epetra::Vector&> (result);
  return applyJacobianTranspose(epetrainput, epetraresult);
}

Abstract::Group::ReturnType Group::applyJacobianTranspose(const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian())
    return Abstract::Group::BadDependency;

  bool status = sharedLinearSystem.getObject()->applyJacobianTranspose(input, result);

  return status == true ? Abstract::Group::Ok : Abstract::Group::Failed;
}

Abstract::Group::ReturnType Group::applyRightPreconditioning(
                      bool useTranspose,
                      Teuchos::ParameterList& params,
                      const Abstract::Vector& input,
                      Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetraInput = dynamic_cast<const NOX::Epetra::Vector&>(input);
  NOX::Epetra::Vector& epetraResult = dynamic_cast<NOX::Epetra::Vector&>(result);

  return applyRightPreconditioning(useTranspose, params, epetraInput, epetraResult);
}

Abstract::Group::ReturnType Group::applyRightPreconditioning(
                       bool useTranspose,
                       Teuchos::ParameterList& linearSolverParams,
                       const NOX::Epetra::Vector& input,
                       NOX::Epetra::Vector& result) const
{

  bool success = false;

  if (!isPreconditioner()) {
    sharedLinearSystem.getObject(this)->destroyPreconditioner();
    sharedLinearSystem.getObject(this)->
      createPreconditioner(xVector, linearSolverParams, false);
    isValidPreconditioner = true;
  }

  success = sharedLinearSystem.getObject()->
    applyRightPreconditioning(useTranspose, linearSolverParams, input, result);

  if (success == true)
    return Abstract::Group::Ok;
  else
    return Abstract::Group::Failed;
}

bool Group::isF() const
{
  return isValidRHS;
}

bool Group::isJacobian() const
{
  return ((sharedLinearSystem.isOwner(this)) && (isValidJacobian));
}

bool Group::isGradient() const
{
  return isValidGrad;
}

bool Group::isNewton() const
{
  return isValidNewton;
}

bool Group::isNormNewtonSolveResidual() const
{
  return isValidNormNewtonSolveResidual;
}

bool Group::isPreconditioner() const
{
  return ((sharedLinearSystem.isOwner(this)) && (isValidPreconditioner) &&
      (sharedLinearSystem.getObject(this)->isPreconditionerConstructed()));
}

bool Group::isConditionNumber() const
{
  return isValidConditionNumber;
}

const Abstract::Vector& Group::getX() const
{
  return xVector;
}

const Abstract::Vector& Group::getF() const
{
  if (!isF()) {
    std::cerr << "ERROR: NOX::Epetra::Group::getF() - invalid RHS" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return RHSVector;
}

double Group::getNormF() const
{
  if (!isF()) {
    std::cerr << "ERROR: NOX::Epetra::Group::getNormF() - invalid RHS" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return RHSVectorPtr->norm();
}

const Abstract::Vector& Group::getGradient() const
{
  if (!isGradient()) {
    std::cerr << "ERROR: NOX::Epetra::Group::getGradient() - invalid gradient" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return gradVector;
}

const Abstract::Vector& Group::getNewton() const
{
  if (!isNewton()) {
    std::cerr << "ERROR: NOX::Epetra::Group::getNewton() - invalid Newton vector" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return NewtonVector;
}

Abstract::Group::ReturnType NOX::Epetra::Group::getNormLastLinearSolveResidual(double& residual) const
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual) {
    residual = normNewtonSolveResidual;
    return NOX::Abstract::Group::Ok;
  }

  // Otherwise give warning since a Newton direction has not been calculated
  // wrt this solution group
  if (utils.isPrintType(Utils::Warning)) {
    std::cout << "ERROR: NOX::Epetra::Group::getNormLastLinearSolveResidual() - "
     << "Group has not performed a Newton solve corresponding to this "
     << "solution vector, or disableLinearSolveResidual(true) was set!" << std::endl;
  }
  return NOX::Abstract::Group::BadDependency;
}

Teuchos::RCP<NOX::Epetra::Interface::Required> Group::
getRequiredInterface()
{
  return userInterfacePtr;
}

Teuchos::RCP<const NOX::Epetra::LinearSystem> Group::
getLinearSystem() const
{
  return sharedLinearSystem.getObject();
}

Teuchos::RCP<NOX::Epetra::LinearSystem> Group::getLinearSystem()
{
  return sharedLinearSystem.getObject(this);
}

void Group::logLastLinearSolveStats(NOX::SolverStats& stats) const
{
  stats.linearSolve.logLinearSolve(linearSolveConverged,
                                   numIterations,
                                   achievedTol,0.0,0.0);
}

bool Group::computeNormNewtonSolveResidual ()
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual)
    return true;

  // Make sure NewtonVector and RHSVector are valid
  // We could return false, but for now we will throw errors
  if (!isValidRHS) {
    std::cerr << "ERROR: NOX::Epetra::Group::computeNormNewtonSolveResidual() - invalid RHS"
     << std::endl;
    throw std::runtime_error("NOX Error");
  }
  if (!isValidNewton) {
    std::cerr << "ERROR: NOX::Epetra::Group::computeNormNewtonSolveResidual() - invalid "
     << "Newton direction" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Allocate the tmpVectorPtr if not already done (deleted in ~Group)
  if (Teuchos::is_null(tmpVectorPtr)) {
    tmpVectorPtr =
      Teuchos::rcp(new Epetra_Vector(RHSVector.getEpetraVector()));
  }
  NOX::Epetra::Vector tmpNoxVector(*tmpVectorPtr, ShapeCopy);

  sharedLinearSystem.getObject()->applyJacobian(NewtonVector, tmpNoxVector);
  tmpNoxVector.update(1.0, RHSVector, 1.0);
  normNewtonSolveResidual = tmpNoxVector.norm();

  isValidNormNewtonSolveResidual = true;

  return true;
}

Abstract::Group::ReturnType NOX::Epetra::Group::
computeJacobianConditionNumber(int maxIters, double tolerance,
                   int krylovSubspaceSize, bool printOutput)
{
  if (!isConditionNumber()) {
    if (!isJacobian()) {
      std::cerr << "ERROR: NOX::Epetra::Group::computeJacobianConditionNumber()"
       << " - Jacobian is invalid wrt the solution." << std::endl;
      throw std::runtime_error("NOX Error");
    }

    if (Teuchos::is_null(azConditionNumberPtr))
      azConditionNumberPtr = Teuchos::rcp(new AztecOOConditionNumber);

    azConditionNumberPtr->
      initialize(*(sharedLinearSystem.getObject()->getJacobianOperator()),
         AztecOOConditionNumber::GMRES_, krylovSubspaceSize,
         printOutput);

    azConditionNumberPtr->computeConditionNumber(maxIters, tolerance);

    conditionNumber = azConditionNumberPtr->getConditionNumber();

    isValidConditionNumber = true;
  }
  return NOX::Abstract::Group::Ok;
}

double NOX::Epetra::Group::getJacobianConditionNumber() const
{
  if (!isConditionNumber()) {
    std::cerr << "ERROR: NOX::Epetra::Group::getJacobianConditionNumber()"
     << " - condition number has not yet been computed!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return conditionNumber;
}

void  NOX::Epetra::Group::disableLinearResidualComputation(const bool disableChoice)
{
  linearResidCompDisabled = disableChoice;
}
