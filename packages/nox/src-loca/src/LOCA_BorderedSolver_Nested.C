// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_Nested.H"
#include "LOCA_BorderedSolver_BorderedOperator.H"
#include "LOCA_BorderedSolver_JacobianOperator.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_BorderedSystem_AbstractGroup.H"
#include "Teuchos_ParameterList.hpp"

LOCA::BorderedSolver::Nested::Nested(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& slvrParams):
  globalData(global_data),
  solverParams(slvrParams),
  solver(),
  grp(),
  unbordered_grp(),
  myWidth(0),
  underlyingWidth(0),
  numConstraints(0)
{
  // Get "Nested Solver" sublist
  Teuchos::RCP<Teuchos::ParameterList> nestedSolverList =
    Teuchos::rcp(&(solverParams->sublist("Nested Bordered Solver")),false);

  // Instantiate underlying solver
  solver =
    globalData->locaFactory->createBorderedSolverStrategy(topParams,
                              nestedSolverList);
}

LOCA::BorderedSolver::Nested::~Nested()
{
}

void
LOCA::BorderedSolver::Nested::setMatrixBlocks(
         const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& oper,
     const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
     const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
     const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  std::string callingFunction =
    "LOCA::BorderedSolver::Nested::setMatrixBlocks()";

  // Cast oper to a bordered operator
  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> op =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::JacobianOperator>(oper);
  if (op == Teuchos::null)
    globalData->locaErrorCheck->throwError(
      callingFunction,
      std::string("Operaror must be of type LOCA::BorderedSolver::JacobianOperator")
      + std::string(" in order to use nested bordered solver strategy."));

  // Get bordered group
  grp = Teuchos::rcp_dynamic_cast<const LOCA::BorderedSystem::AbstractGroup>(op->getGroup());
  if (grp == Teuchos::null)
    globalData->locaErrorCheck->throwError(
      callingFunction,
      std::string("Group must be of type LOCA::BorderedSystem::AbstractGroup")
      + std::string(" in order to use nested bordered solver strategy."));

  Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> con_mvdx = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (con_mvdx == Teuchos::null)
    globalData->locaErrorCheck->throwError(
         callingFunction,
         "Constraints object must be of type ConstraintInterfaceMVDX");

  bool isZeroA = (blockA.get() == NULL);
  bool isZeroB = con_mvdx->isDXZero();
  bool isZeroC = (blockC.get() == NULL);
  Teuchos::RCP<const NOX::Abstract::MultiVector> blockB_dx;
  if (!isZeroB)
    blockB_dx = Teuchos::rcp(con_mvdx->getDX(), false);

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

  // Get unbordered group
  unbordered_grp = grp->getUnborderedGroup();

  // get number of outer constraints
  if (isZeroB)
    numConstraints = blockC->numRows();
  else
    numConstraints = blockB_dx->numVectors();

  // Get total bordered width
  underlyingWidth = grp->getBorderedWidth();
  myWidth = underlyingWidth + numConstraints;

  // combine blocks
  bool isCombinedAZero = grp->isCombinedAZero();
  bool isCombinedBZero = grp->isCombinedBZero();
  bool isCombinedCZero = grp->isCombinedCZero();
  Teuchos::RCP<NOX::Abstract::MultiVector> A;
  Teuchos::RCP<NOX::Abstract::MultiVector> B;
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> C;

  if (!isCombinedAZero || !isZeroA) {
    A = unbordered_grp->getX().createMultiVector(myWidth);
    A->init(0.0);
  }
  if (!isCombinedBZero || !isZeroB) {
    B = unbordered_grp->getX().createMultiVector(myWidth);
    B->init(0.0);

  }
  if (!isCombinedCZero || !isZeroC) {
    C = Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(myWidth,
                                 myWidth));
    C->putScalar(0.0);
  }

  std::vector<int> idx1(underlyingWidth);
  for (int i=0; i<underlyingWidth; i++)
    idx1[i] = i;
  if (!isCombinedAZero) {
    Teuchos::RCP<NOX::Abstract::MultiVector> underlyingA =
      A->subView(idx1);
    grp->fillA(*underlyingA);
  }
  if (!isCombinedBZero) {
    Teuchos::RCP<NOX::Abstract::MultiVector> underlyingB =
      B->subView(idx1);
    grp->fillB(*underlyingB);
  }
  if (!isCombinedCZero) {
    NOX::Abstract::MultiVector::DenseMatrix underlyingC(Teuchos::View,
                            *C,
                            underlyingWidth,
                            underlyingWidth,
                            0, 0);
    grp->fillC(underlyingC);
  }

  std::vector<int> idx2(numConstraints);
  for (int i=0; i<numConstraints; i++)
    idx2[i] = underlyingWidth+i;
  if (!isZeroA) {
    Teuchos::RCP<NOX::Abstract::MultiVector> my_A_x = A->subView(idx2);
    NOX::Abstract::MultiVector::DenseMatrix my_A_p(Teuchos::View, *C,
                           underlyingWidth,
                           numConstraints, 0,
                           underlyingWidth);
    grp->extractSolutionComponent(*blockA, *my_A_x);
    grp->extractParameterComponent(false, *blockA, my_A_p);
  }

  if (!isZeroB) {
    Teuchos::RCP<NOX::Abstract::MultiVector> my_B_x = B->subView(idx2);
    NOX::Abstract::MultiVector::DenseMatrix my_B_p(Teuchos::View, *C,
                           numConstraints,
                           underlyingWidth,
                           underlyingWidth, 0);
    grp->extractSolutionComponent(*blockB_dx, *my_B_x);
    grp->extractParameterComponent(true, *blockB_dx, my_B_p);
  }

  if (!isZeroC) {
    NOX::Abstract::MultiVector::DenseMatrix my_CC(Teuchos::View, *C,
                          numConstraints,
                          numConstraints,
                          underlyingWidth,
                          underlyingWidth);
    my_CC.assign(*blockC);
  }

  // Create unbordered operator
  Teuchos::RCP<LOCA::BorderedSolver::AbstractOperator> unbordered_op =
    Teuchos::rcp(new LOCA::BorderedSolver::JacobianOperator(unbordered_grp));

  // set blocks in solver
  solver->setMatrixBlocksMultiVecConstraint(unbordered_op, A, B, C);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::initForSolve()
{
  return solver->initForSolve();
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::initForTransposeSolve()
{
  return solver->initForTransposeSolve();
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::apply(
              const NOX::Abstract::MultiVector& X,
              const NOX::Abstract::MultiVector::DenseMatrix& Y,
              NOX::Abstract::MultiVector& U,
              NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  int num_cols = X.numVectors();
  Teuchos::RCP<NOX::Abstract::MultiVector> XX =
    unbordered_grp->getX().createMultiVector(num_cols);
  Teuchos::RCP<NOX::Abstract::MultiVector> UU =
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix VV(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
                          underlyingWidth, num_cols,
                          0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
                          numConstraints, num_cols,
                          underlyingWidth, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV1(Teuchos::View, VV,
                          underlyingWidth, num_cols,
                          0, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV2(Teuchos::View, VV,
                          numConstraints, num_cols,
                          underlyingWidth, 0);

  grp->extractSolutionComponent(X, *XX);
  grp->extractParameterComponent(false, X, YY1);
  YY2.assign(Y);

  NOX::Abstract::Group::ReturnType status =
    solver->apply(*XX, YY, *UU, VV);

  V.assign(VV2);
  grp->loadNestedComponents(*UU, VV1, U);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::applyTranspose(
              const NOX::Abstract::MultiVector& X,
              const NOX::Abstract::MultiVector::DenseMatrix& Y,
              NOX::Abstract::MultiVector& U,
              NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  int num_cols = X.numVectors();
  Teuchos::RCP<NOX::Abstract::MultiVector> XX =
    unbordered_grp->getX().createMultiVector(num_cols);
  Teuchos::RCP<NOX::Abstract::MultiVector> UU =
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix VV(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
                          underlyingWidth, num_cols,
                          0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
                          numConstraints, num_cols,
                          underlyingWidth, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV1(Teuchos::View, VV,
                          underlyingWidth, num_cols,
                          0, 0);
  NOX::Abstract::MultiVector::DenseMatrix VV2(Teuchos::View, VV,
                          numConstraints, num_cols,
                          underlyingWidth, 0);

  grp->extractSolutionComponent(X, *XX);
  grp->extractParameterComponent(false, X, YY1);
  YY2.assign(Y);

  NOX::Abstract::Group::ReturnType status =
    solver->applyTranspose(*XX, YY, *UU, VV);

  V.assign(VV2);
  grp->loadNestedComponents(*UU, VV1, U);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::applyInverse(
                  Teuchos::ParameterList& params,
                  const NOX::Abstract::MultiVector* F,
                  const NOX::Abstract::MultiVector::DenseMatrix* G,
                  NOX::Abstract::MultiVector& X,
                  NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
  }

  int num_cols = X.numVectors();
  Teuchos::RCP<NOX::Abstract::MultiVector> FF;
  if (!isZeroF)
    FF = unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix GG(myWidth, num_cols);
  GG.putScalar(0.0);

  if (!isZeroF) {
    NOX::Abstract::MultiVector::DenseMatrix GG1(Teuchos::View, GG,
                        underlyingWidth, num_cols,
                        0, 0);
    grp->extractSolutionComponent(*F, *FF);
    grp->extractParameterComponent(false, *F, GG1);
  }
  if (!isZeroG) {
    NOX::Abstract::MultiVector::DenseMatrix GG2(Teuchos::View, GG,
                        numConstraints, num_cols,
                        underlyingWidth, 0);
    GG2.assign(*G);
  }

  Teuchos::RCP<NOX::Abstract::MultiVector> XX =
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
                          underlyingWidth, num_cols,
                          0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
                          numConstraints, num_cols,
                          underlyingWidth, 0);

  NOX::Abstract::Group::ReturnType status =
    solver->applyInverse(params, FF.get(), &GG, *XX, YY);

  Y.assign(YY2);
  grp->loadNestedComponents(*XX, YY1, X);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::applyInverseTranspose(
                  Teuchos::ParameterList& params,
                  const NOX::Abstract::MultiVector* F,
                  const NOX::Abstract::MultiVector::DenseMatrix* G,
                  NOX::Abstract::MultiVector& X,
                  NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
  }

  int num_cols = X.numVectors();
  Teuchos::RCP<NOX::Abstract::MultiVector> FF;
  if (!isZeroF)
    FF = unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix GG(myWidth, num_cols);
  GG.putScalar(0.0);

  if (!isZeroF) {
    NOX::Abstract::MultiVector::DenseMatrix GG1(Teuchos::View, GG,
                        underlyingWidth, num_cols,
                        0, 0);
    grp->extractSolutionComponent(*F, *FF);
    grp->extractParameterComponent(false, *F, GG1);
  }
  if (!isZeroG) {
    NOX::Abstract::MultiVector::DenseMatrix GG2(Teuchos::View, GG,
                        numConstraints, num_cols,
                        underlyingWidth, 0);
    GG2.assign(*G);
  }

  Teuchos::RCP<NOX::Abstract::MultiVector> XX =
    unbordered_grp->getX().createMultiVector(num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY(myWidth, num_cols);
  NOX::Abstract::MultiVector::DenseMatrix YY1(Teuchos::View, YY,
                          underlyingWidth, num_cols,
                          0, 0);
  NOX::Abstract::MultiVector::DenseMatrix YY2(Teuchos::View, YY,
                          numConstraints, num_cols,
                          underlyingWidth, 0);

  NOX::Abstract::Group::ReturnType status =
    solver->applyInverseTranspose(params, FF.get(), &GG, *XX, YY);

  Y.assign(YY2);
  grp->loadNestedComponents(*XX, YY1, X);

  return status;
}

