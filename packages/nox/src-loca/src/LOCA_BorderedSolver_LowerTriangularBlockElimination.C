// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_LowerTriangularBlockElimination.H"
#include "LOCA_BorderedSolver_AbstractOperator.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "LOCA_MultiContinuation_MultiVecConstraint.H"
#include "Teuchos_LAPACK.hpp"    // for LAPACK solve

LOCA::BorderedSolver::LowerTriangularBlockElimination::
LowerTriangularBlockElimination(
     const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::BorderedSolver::LowerTriangularBlockElimination::
~LowerTriangularBlockElimination()
{
}


NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::
solve(Teuchos::ParameterList& params,
      const LOCA::BorderedSolver::AbstractOperator& op,
      const LOCA::MultiContinuation::ConstraintInterface& B,
      const NOX::Abstract::MultiVector::DenseMatrix& C,
      const NOX::Abstract::MultiVector* F,
      const NOX::Abstract::MultiVector::DenseMatrix* G,
      NOX::Abstract::MultiVector& X,
      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  std::string callingFunction =
    "LOCA::BorderedSolver::LowerTriangularBlockElimination::solve()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Determine if X or Y is zero
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);
  bool isZeroB = B.isDXZero();
  bool isZeroX = isZeroF;
  bool isZeroY = isZeroG && (isZeroB  || isZeroX);

  // First compute X
  if (isZeroX)
    X.init(0.0);
  else {
    // Solve X = J^-1 F, note F must be nonzero
    status = op.applyInverse(params, *F, X);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Now compute Y
  if (isZeroY)
    Y.putScalar(0.0);
  else {
    // Compute G - B^T*X and store in Y
    if (isZeroG)
      B.multiplyDX(-1.0, X, Y);
    else {
      Y.assign(*G);
      if (!isZeroB && !isZeroX) {
    NOX::Abstract::MultiVector::DenseMatrix T(Y.numRows(),Y.numCols());
    B.multiplyDX(1.0, X, T);
    Y -= T;
      }
    }

    // Overwrite Y with Y = C^-1 * (G - B^T*X)
    NOX::Abstract::MultiVector::DenseMatrix M(Teuchos::Copy, C);
    int *ipiv = new int[M.numRows()];
    Teuchos::LAPACK<int,double> L;
    int info;
    L.GETRF(M.numRows(), M.numCols(), M.values(), M.stride(), ipiv, &info);
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(
                                  status,
                                  finalStatus,
                                  callingFunction);
    }
    L.GETRS('N', M.numRows(), Y.numCols(), M.values(), M.stride(), ipiv,
        Y.values(), Y.stride(), &info);
    delete [] ipiv;
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(
                                 status,
                                 finalStatus,
                                 callingFunction);
    }
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::
solve(Teuchos::ParameterList& params,
      const LOCA::BorderedSolver::AbstractOperator& op,
      const NOX::Abstract::MultiVector& B,
      const NOX::Abstract::MultiVector::DenseMatrix& C,
      const NOX::Abstract::MultiVector* F,
      const NOX::Abstract::MultiVector::DenseMatrix* G,
      NOX::Abstract::MultiVector& X,
      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  // Create a constraint out of B
  LOCA::MultiContinuation::MultiVecConstraint cB(Teuchos::rcp(&B,false));

  // Call other version
  return solve(params, op, cB, C, F, G, X, Y);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::
solveTranspose(Teuchos::ParameterList& params,
           const LOCA::BorderedSolver::AbstractOperator& op,
           const LOCA::MultiContinuation::ConstraintInterface& B,
           const NOX::Abstract::MultiVector::DenseMatrix& C,
           const NOX::Abstract::MultiVector* F,
           const NOX::Abstract::MultiVector::DenseMatrix* G,
           NOX::Abstract::MultiVector& X,
           NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  std::string callingFunction =
    "LOCA::BorderedSolver::LowerTriangularBlockElimination::solveTranspose()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Determine if X or Y is zero
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);
  bool isZeroB = B.isDXZero();
  bool isZeroX = isZeroF;
  bool isZeroY = isZeroG && (isZeroB  || isZeroX);

  // First compute X
  if (isZeroX)
    X.init(0.0);
  else {
    // Solve X = J^-T F, note F must be nonzero
    status = op.applyInverseTranspose(params, *F, X);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Now compute Y
  if (isZeroY)
    Y.putScalar(0.0);
  else {
    // Compute G - B^T*X and store in Y
    if (isZeroG)
      B.multiplyDX(-1.0, X, Y);
    else {
      Y.assign(*G);
      if (!isZeroB && !isZeroX) {
    NOX::Abstract::MultiVector::DenseMatrix T(Y.numRows(),Y.numCols());
    B.multiplyDX(1.0, X, T);
    Y -= T;
      }
    }

    // Overwrite Y with Y = C^-T * (G - B^T*X)
    NOX::Abstract::MultiVector::DenseMatrix M(Teuchos::Copy, C);
    int *ipiv = new int[M.numRows()];
    Teuchos::LAPACK<int,double> L;
    int info;
    L.GETRF(M.numRows(), M.numCols(), M.values(), M.stride(), ipiv, &info);
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(
                                  status,
                                  finalStatus,
                                  callingFunction);
    }
    L.GETRS('T', M.numRows(), Y.numCols(), M.values(), M.stride(), ipiv,
        Y.values(), Y.stride(), &info);
    delete [] ipiv;
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(
                                 status,
                                 finalStatus,
                                 callingFunction);
    }
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::
solveTranspose(Teuchos::ParameterList& params,
           const LOCA::BorderedSolver::AbstractOperator& op,
           const NOX::Abstract::MultiVector& B,
           const NOX::Abstract::MultiVector::DenseMatrix& C,
           const NOX::Abstract::MultiVector* F,
           const NOX::Abstract::MultiVector::DenseMatrix* G,
           NOX::Abstract::MultiVector& X,
           NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  // Create a constraint out of B
  LOCA::MultiContinuation::MultiVecConstraint cB(Teuchos::rcp(&B,false));

  // Call other version
  return solveTranspose(params, op, cB, C, F, G, X, Y);
}
