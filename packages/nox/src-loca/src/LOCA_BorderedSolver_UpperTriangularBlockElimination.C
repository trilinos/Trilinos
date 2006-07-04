// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "LOCA_BorderedSolver_UpperTriangularBlockElimination.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_LAPACK.hpp"    // for LAPACK solve
#include "LOCA_Abstract_TransposeSolveGroup.H"

LOCA::BorderedSolver::UpperTriangularBlockElimination::
UpperTriangularBlockElimination(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::BorderedSolver::UpperTriangularBlockElimination::
~UpperTriangularBlockElimination()
{
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::UpperTriangularBlockElimination::
solve(Teuchos::ParameterList& params,
      const NOX::Abstract::Group& grp,
      const NOX::Abstract::MultiVector* A,
      const NOX::Abstract::MultiVector::DenseMatrix& C,
      const NOX::Abstract::MultiVector* F,
      const NOX::Abstract::MultiVector::DenseMatrix* G,
      NOX::Abstract::MultiVector& X,
      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
 string callingFunction = 
    "LOCA::BorderedSolver::UpperTriangularBlockElimination::solve()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Determine if X or Y is zero
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);
  bool isZeroA = (A == NULL);
  bool isZeroY = isZeroG;
  bool isZeroX = isZeroF && (isZeroA || isZeroY);

  // First compute Y
  if (isZeroY)
    Y.putScalar(0.0);
  else {
    // Solve Y = C^-1 * G
    NOX::Abstract::MultiVector::DenseMatrix M(C);
    int *ipiv = new int[M.numRows()];
    Teuchos::LAPACK<int,double> L;
    int info;
    
    Y.assign(*G);
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

  // Now compute X
  if (isZeroX)
    X.init(0.0);
  else if (isZeroA || isZeroY) {
    // Solve X = J^-1 F, note F must be nonzero
    status = grp.applyJacobianInverseMultiVector(params, *F, X);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  else {
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS;

    if (isZeroF) {
      RHS = A->clone(Y.numCols());
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 0.0);
    }
    else {
      RHS = F->clone(NOX::DeepCopy);
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 1.0);
    }
    // Solve X = J^-1 (F-A*Y)
    status = grp.applyJacobianInverseMultiVector(params, *RHS, X);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::UpperTriangularBlockElimination::
solveTranspose(Teuchos::ParameterList& params,
	       const LOCA::Abstract::TransposeSolveGroup& grp,
	       const NOX::Abstract::MultiVector* A,
	       const NOX::Abstract::MultiVector::DenseMatrix& C,
	       const NOX::Abstract::MultiVector* F,
	       const NOX::Abstract::MultiVector::DenseMatrix* G,
	       NOX::Abstract::MultiVector& X,
	       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
 string callingFunction = 
    "LOCA::BorderedSolver::UpperTriangularBlockElimination::solveTranspose()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Determine if X or Y is zero
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);
  bool isZeroA = (A == NULL);
  bool isZeroY = isZeroG;
  bool isZeroX = isZeroF && (isZeroA || isZeroY);

  // First compute Y
  if (isZeroY)
    Y.putScalar(0.0);
  else {
    // Solve Y = C^-T * G
    NOX::Abstract::MultiVector::DenseMatrix M(C);
    int *ipiv = new int[M.numRows()];
    Teuchos::LAPACK<int,double> L;
    int info;
    
    Y.assign(*G);
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

  // Now compute X
  if (isZeroX)
    X.init(0.0);
  else if (isZeroA || isZeroY) {
    // Solve X = J^-T F, note F must be nonzero
    status = grp.applyJacobianTransposeInverseMultiVector(params, *F, X);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  else {
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS;

    if (isZeroF) {
      RHS = A->clone(Y.numCols());
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 0.0);
    }
    else {
      RHS = F->clone(NOX::DeepCopy);
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 1.0);
    }
    // Solve X = J^-T (F-A*Y)
    status = grp.applyJacobianTransposeInverseMultiVector(params, *RHS, X);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}
