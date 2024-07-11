// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_EpetraAugmented.H"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "NOX_Epetra_MultiVector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_Epetra_Group.H"
#include "LOCA_Epetra_AugmentedOp.H"
#include "LOCA_BorderedSolver_LowerTriangularBlockElimination.H"
#include "LOCA_BorderedSolver_UpperTriangularBlockElimination.H"
#include "LOCA_BorderedSolver_JacobianOperator.H"

LOCA::BorderedSolver::EpetraAugmented::EpetraAugmented(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
     const Teuchos::RCP<Teuchos::ParameterList>& slvrParams):
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  op(),
  A(),
  B(),
  C(),
  constraints(),
  numConstraints(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true)
{
}

LOCA::BorderedSolver::EpetraAugmented::~EpetraAugmented()
{
}

void
LOCA::BorderedSolver::EpetraAugmented::setMatrixBlocks(
         const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op_,
     const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
     const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
     const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  std::string callingFunction =
    "LOCA::BorderedSolver::EpetraAugmented::setMatrixBlocks";

  op = op_;

  // Get Jacobian operator
  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> jacOp =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::JacobianOperator>(op);

  Teuchos::RCP<const NOX::Abstract::Group> group =
    jacOp->getGroup();

  // Cast away const
  Teuchos::RCP<NOX::Abstract::Group> non_const_group =
    Teuchos::rcp_const_cast<NOX::Abstract::Group>(group);

  // Cast group to an Epetra group
  grp = Teuchos::rcp_dynamic_cast<LOCA::Epetra::Group>(non_const_group);
  if (grp.get() == NULL)
    globalData->locaErrorCheck->throwError(
                    callingFunction,
                    "Group object must be an Epetra group");

  A = blockA;

  // Cast constraints to a ConstraintInterfaceMVDX
  constraints = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (constraints.get() == NULL)
    globalData->locaErrorCheck->throwError(
         callingFunction,
         "Constraints object must be of type ConstraintInterfaceMVDX");
  C = blockC;

  // Determine which blocks are zero
  isZeroA = (A.get() == NULL);
  isZeroB = constraints->isDXZero();
  isZeroC = (C.get() == NULL);

  // Get multivec constraint if it is nonzero
  if (isZeroB)
    B = Teuchos::null;
  else
    B = Teuchos::rcp(constraints->getDX(), false);

  // ensure blocks B and C are not both zero
  if (isZeroB && isZeroC)
    globalData->locaErrorCheck->throwError(
                        callingFunction,
                        "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC)
    globalData->locaErrorCheck->throwError(
                         callingFunction,
                         "Blocks A and C cannot both be zero");

  if (isZeroB)
    numConstraints = C->numRows();
  else
    numConstraints = B->numVectors();

  // We only use the augmented technique if A and B are nonzero
  // (otherwise we use a block elimination).  If C is zero, we just create
  // a matrix of zeros and apply the standard algorithm
  if (isZeroC && !isZeroA && !isZeroB) {

    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tmpC =
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
                               B->numVectors(),
                               B->numVectors()));
    tmpC->putScalar(0.0);
    C = tmpC;
    isZeroC = false;
  }

}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::initForSolve()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::initForTransposeSolve()
{
  return NOX::Abstract::Group::Ok;
}


NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::apply(
              const NOX::Abstract::MultiVector& X,
              const NOX::Abstract::MultiVector::DenseMatrix& Y,
              NOX::Abstract::MultiVector& U,
              NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  // Compute J*X
  NOX::Abstract::Group::ReturnType status =
    grp->applyJacobianMultiVector(X, U);

  // Compute J*X + A*Y
  if (!isZeroA)
    U.update(Teuchos::NO_TRANS, 1.0, *A, Y, 1.0);

  // Compute B^T*X
  if (!isZeroB)
    constraints->multiplyDX(1.0, X, V);

  // Compute B^T*X + C*Y
  if (!isZeroC) {
    int e;
    if (isZeroB)
      e = V.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 0.0);
    else
      e = V.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 1.0);
    if (e < 0)
      status = NOX::Abstract::Group::Failed;
  }

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::applyTranspose(
              const NOX::Abstract::MultiVector& X,
              const NOX::Abstract::MultiVector::DenseMatrix& Y,
              NOX::Abstract::MultiVector& U,
              NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  // Compute J*X
  NOX::Abstract::Group::ReturnType status =
    grp->applyJacobianTransposeMultiVector(X, U);

  // Compute J*X + B*Y
  if (!isZeroA)
    constraints->addDX(Teuchos::NO_TRANS, 1.0, Y, 1.0, U);

  // Compute A^T*X
  if (!isZeroB)
    X.multiply(1.0, *A, V);

  // Compute A^T*X + C^T*Y
  if (!isZeroC) {
    int e;
    if (isZeroB)
      e = V.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 0.0);
    else
      e = V.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 1.0);
    if (e < 0)
      status = NOX::Abstract::Group::Failed;
  }

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::applyInverse(
                  Teuchos::ParameterList& params,
                  const NOX::Abstract::MultiVector* F,
                  const NOX::Abstract::MultiVector::DenseMatrix* G,
                  NOX::Abstract::MultiVector& X,
                  NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  std::string callingFunction =
    "LOCA::BorderedSolver::EpetraAugmented::applyInverse()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F and G are zero, return zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  if (isZeroA) {
    LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
    return ltbe.solve(params, *op, *constraints, *C, F, G, X, Y);
  }

  if (isZeroB) {
    LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
    return utbe.solve(params, *op, A.get(), *C, F, G, X, Y);
  }

  // Get underlying Epetra vectors
  Teuchos::RCP<const NOX::Epetra::MultiVector> nox_epetra_a =
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(A);
  Teuchos::RCP<const NOX::Epetra::MultiVector> nox_epetra_b =
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(B);
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_x =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(Teuchos::rcp(&X,
                                     false));
  Teuchos::RCP<const NOX::Epetra::MultiVector> nox_epetra_f =
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(Teuchos::rcp(F, false));
  Teuchos::RCP<const Epetra_MultiVector> epetra_a =
    Teuchos::rcp(&(nox_epetra_a->getEpetraMultiVector()), false);
  Teuchos::RCP<const Epetra_MultiVector> epetra_b =
    Teuchos::rcp(&(nox_epetra_b->getEpetraMultiVector()), false);
  Teuchos::RCP<Epetra_MultiVector> epetra_x =
    Teuchos::rcp(&(nox_epetra_x->getEpetraMultiVector()), false);
  Teuchos::RCP<const Epetra_MultiVector> epetra_f =
     Teuchos::rcp(&(nox_epetra_f->getEpetraMultiVector()), false);

  const NOX::Epetra::Vector& solution_vec =
    dynamic_cast<const NOX::Epetra::Vector&>(grp->getX());

  // Get linear system
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
    grp->getLinearSystem();

  // Get Jacobian
  Teuchos::RCP<Epetra_Operator> jac =
    linSys->getJacobianOperator();

  // Set Jacobian
  linSys->setJacobianOperatorForSolve(jac);

  // Create the preconditioner
  linSys->destroyPreconditioner();
  linSys->createPreconditioner(solution_vec, params, false);

  // Get preconditioner
  Teuchos::RCP<Epetra_Operator> prec =
    linSys->getGeneratedPrecOperator();

  // Create augmented operators
  LOCA::Epetra::AugmentedOp extended_jac(globalData, jac, epetra_a,
                     epetra_b, C);
  LOCA::Epetra::AugmentedOp extended_prec(globalData, prec, epetra_a,
                      epetra_b, C);

  // Set augmented operator in linear system
  linSys->setJacobianOperatorForSolve(Teuchos::rcp(&extended_jac,false));
  linSys->setPrecOperatorForSolve(Teuchos::rcp(&extended_prec,false));

  // Create augmented Epetra vectors for x, f
  Teuchos::RCP<Epetra_MultiVector> epetra_augmented_x =
    extended_jac.buildEpetraAugmentedMultiVec(*epetra_x, &Y, false);
  Teuchos::RCP<Epetra_MultiVector> epetra_augmented_f =
    extended_jac.buildEpetraAugmentedMultiVec(*epetra_f, G, true);

  // Create augmented NOX::Epetra::MultiVectors as views
  NOX::Epetra::MultiVector nox_epetra_augmented_x(
                     epetra_augmented_x,
                     NOX::DeepCopy,
                     NOX::Epetra::MultiVector::CreateView);
  NOX::Epetra::MultiVector nox_epetra_augmented_f(
                     epetra_augmented_f,
                     NOX::DeepCopy,
                     NOX::Epetra::MultiVector::CreateView);

  // Solve for each RHS
  int m = nox_epetra_augmented_x.numVectors();
  for (int i=0; i<m; i++) {
    extended_jac.init(*((*epetra_f)(i)));
    extended_prec.init(*((*epetra_f)(i)));
    bool stat =
      linSys->applyJacobianInverse(
        params,
        dynamic_cast<NOX::Epetra::Vector&>(nox_epetra_augmented_f[i]),
        dynamic_cast<NOX::Epetra::Vector&>(nox_epetra_augmented_x[i]));
    if (stat == true)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::NotConverged;
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Set results
  extended_jac.init(*epetra_x);
  extended_jac.setEpetraAugmentedMultiVec(*epetra_x, Y,
                      *epetra_augmented_x);

  // Set original Jacobian in linear system
  linSys->setJacobianOperatorForSolve(jac);
  linSys->destroyPreconditioner();
  //linSys.setPrecOperatorForSolve(*prec);

  return finalStatus;
}
