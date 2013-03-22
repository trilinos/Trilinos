// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_BorderedSolver_UpperTriangularBlockElimination.H"
#include "LOCA_BorderedSolver_AbstractOperator.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_LAPACK.hpp"    // for LAPACK solve
#include "LOCA_Abstract_TransposeSolveGroup.H"

LOCA::BorderedSolver::UpperTriangularBlockElimination::
UpperTriangularBlockElimination(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data) : 
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
      const LOCA::BorderedSolver::AbstractOperator& op,
      const NOX::Abstract::MultiVector* A,
      const NOX::Abstract::MultiVector::DenseMatrix& C,
      const NOX::Abstract::MultiVector* F,
      const NOX::Abstract::MultiVector::DenseMatrix* G,
      NOX::Abstract::MultiVector& X,
      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
 std::string callingFunction = 
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
    status = op.applyInverse(params, *F, X);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  else {
    Teuchos::RCP<NOX::Abstract::MultiVector> RHS;

    if (isZeroF) {
      RHS = A->clone(Y.numCols());
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 0.0);
    }
    else {
      RHS = F->clone(NOX::DeepCopy);
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 1.0);
    }
    // Solve X = J^-1 (F-A*Y)
    status = op.applyInverse(params, *RHS, X);
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
	       const LOCA::BorderedSolver::AbstractOperator& op,
	       const NOX::Abstract::MultiVector* A,
	       const NOX::Abstract::MultiVector::DenseMatrix& C,
	       const NOX::Abstract::MultiVector* F,
	       const NOX::Abstract::MultiVector::DenseMatrix* G,
	       NOX::Abstract::MultiVector& X,
	       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
 std::string callingFunction = 
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
    status = op.applyInverseTranspose(params, *F, X);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  else {
    Teuchos::RCP<NOX::Abstract::MultiVector> RHS;

    if (isZeroF) {
      RHS = A->clone(Y.numCols());
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 0.0);
    }
    else {
      RHS = F->clone(NOX::DeepCopy);
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 1.0);
    }
    // Solve X = J^-T (F-A*Y)
    status = op.applyInverseTranspose(params, *RHS, X);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}
