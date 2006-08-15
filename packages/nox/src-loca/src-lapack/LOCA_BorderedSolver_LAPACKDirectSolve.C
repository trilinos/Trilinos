// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

#include "LOCA_BorderedSolver_LAPACKDirectSolve.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_LAPACK.H"
#include "Teuchos_LAPACK.hpp"

LOCA::BorderedSolver::LAPACKDirectSolve::LAPACKDirectSolve(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  A(),
  B(),
  C(),
  augmentedSolver(),
  n(0),
  m(0),
  N(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isZeroF(true),
  isZeroG(true)
{
}

LOCA::BorderedSolver::LAPACKDirectSolve::~LAPACKDirectSolve()
{
}

void
LOCA::BorderedSolver::LAPACKDirectSolve::setMatrixBlocks(
         const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  string callingFunction = 
    "LOCA::BorderedSolver::LAPACKDirectSolve::setMatrixBlocks()";
  const NOX::LAPACK::Vector *v;
  const NOX::Abstract::MultiVector *BV;

  // Set block pointers
  grp = Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(group);
  if (grp.get() == NULL)
    globalData->locaErrorCheck->throwError(
	      callingFunction,
	      string("Group argument is not of type LOCA::LAPACK::Group!\n") + 
	      string("The LAPACK Direct Solve bordered solver method can\n") +
	      string("only be used with LAPACK groups."));
  A = blockA;
  
  B = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (B.get() == NULL)
    globalData->locaErrorCheck->throwError(
	     callingFunction,
	     string("Constraint argument is not of type\n") + 
	     string("LOCA::MultiContinuation::ConstraintInterfaceMVDX!\n") +
	     string("The LAPACK Direct Solve bordered solver method can\n") +
	     string("only be used with constraints that support obtaining\n") +
	     string("the constraint derivative as a multivector."));

  C = blockC;

  isZeroA = (A.get() == NULL);
  isZeroB = B->isDXZero();
  isZeroC = (C.get() == NULL);

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

  // Get the Jacobian matrix and size
  const NOX::LAPACK::Matrix<double>& J = grp->getJacobianMatrix();
  n = J.numRows();

  // Get the number of additional rows/columns
  if (!isZeroA)
    m = A->numVectors();
  else
    m = C->numCols();

  // Form a new (n+m) x (n+m) matrix if this is a new size
  if (n+m != N) {
    N = n+m;
    augmentedSolver = Teuchos::rcp(new NOX::LAPACK::LinearSolver<double>(N));
  }
  else {
    augmentedSolver->reset();
  }
  NOX::LAPACK::Matrix<double>& augmentedJ = augmentedSolver->getMatrix();

  // Copy Jacobian
  for (int j=0; j<n; j++)
    for (int i=0; i<n; i++)
      augmentedJ(i,j) = J(i,j);

  // Copy A
  if (isZeroA) {
    for (int j=0; j<m; j++) 
      for (int i=0; i<n; i++)
	augmentedJ(i,j+n) = 0.0;
  }
  else {
    for (int j=0; j<m; j++) {
      v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*A)[j]);
      for (int i=0; i<n; i++)
	augmentedJ(i,j+n) = (*v)(i);
    }
  }

  // Copy B
  if (isZeroB) {
    for (int i=0; i<m; i++) 
      for (int j=0; j<n; j++)
	augmentedJ(i+n,j) = 0.0;
  }
  else {
    BV = B->getDX();
    for (int i=0; i<m; i++) {
      v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*BV)[i]);
      for (int j=0; j<n; j++)
	augmentedJ(i+n,j) = (*v)(j);
    }
  }

  // Copy C
  if (isZeroC) {
    for (int j=0; j<m; j++)
      for (int i=0; i<m; i++)
	augmentedJ(i+n,j+n) = 0.0;
  }
  else {
    for (int j=0; j<m; j++)
      for (int i=0; i<m; i++)
	augmentedJ(i+n,j+n) = (*C)(i,j);
  }
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::LAPACKDirectSolve::initForSolve()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::LAPACKDirectSolve::initForTransposeSolve()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::LAPACKDirectSolve::apply(
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
    B->multiplyDX(1.0, X, V);

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
LOCA::BorderedSolver::LAPACKDirectSolve::applyTranspose(
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
    B->addDX(Teuchos::NO_TRANS, 1.0, Y, 1.0, U);

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
LOCA::BorderedSolver::LAPACKDirectSolve::applyInverse(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F & G are zero, the solution is zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  int numColsRHS;
  const NOX::LAPACK::Vector *v;
  NOX::LAPACK::Vector *w;

  if (!isZeroF)
    numColsRHS = F->numVectors();
  else 
    numColsRHS = G->numCols();
    
  // Concatenate F & G into a single matrix
  NOX::LAPACK::Matrix<double> RHS(N,numColsRHS);
  if (isZeroF) {
    for (int j=0; j<numColsRHS; j++)
      for (int i=0; i<n; i++)
	RHS(i,j) = 0.0;
  }
  else {
    for (int j=0; j<numColsRHS; j++) {
      v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*F)[j]);
      for (int i=0; i<n; i++)
	RHS(i,j) = (*v)(i);
    }
  }
  if (isZeroG) {
    for (int j=0; j<numColsRHS; j++)
      for (int i=0; i<m; i++)
	RHS(i+n,j) = 0.0;
  }
  else {
    for (int j=0; j<numColsRHS; j++)
      for (int i=0; i<m; i++)
	RHS(i+n,j) = (*G)(i,j);
  }

  // Solve for LHS
  bool res = augmentedSolver->solve(false, numColsRHS, &RHS(0,0));

  // Copy result into X and Y
  for (int j=0; j<numColsRHS; j++) {
    w = dynamic_cast<NOX::LAPACK::Vector*>(&X[j]);
    for (int i=0; i<n; i++)
      (*w)(i) = RHS(i,j);
    for (int i=0; i<m; i++)
      Y(i,j) = RHS(n+i,j);
  }

  if (!res)
    return NOX::Abstract::Group::Failed;
  else
    return NOX::Abstract::Group::Ok;
  
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::LAPACKDirectSolve::applyInverseTranspose(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F & G are zero, the solution is zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  int numColsRHS;
  const NOX::LAPACK::Vector *v;
  NOX::LAPACK::Vector *w;

  if (!isZeroF)
    numColsRHS = F->numVectors();
  else 
    numColsRHS = G->numCols();
    
  // Concatenate F & G into a single matrix
  NOX::LAPACK::Matrix<double> RHS(N,numColsRHS);
  if (isZeroF) {
    for (int j=0; j<numColsRHS; j++)
      for (int i=0; i<n; i++)
	RHS(i,j) = 0.0;
  }
  else {
    for (int j=0; j<numColsRHS; j++) {
      v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*F)[j]);
      for (int i=0; i<n; i++)
	RHS(i,j) = (*v)(i);
    }
  }
  if (isZeroG) {
    for (int j=0; j<numColsRHS; j++)
      for (int i=0; i<m; i++)
	RHS(i+n,j) = 0.0;
  }
  else {
    for (int j=0; j<numColsRHS; j++)
      for (int i=0; i<m; i++)
	RHS(i+n,j) = (*G)(i,j);
  }

  // Solve for LHS
  bool res = augmentedSolver->solve(true, numColsRHS, &RHS(0,0));

  // Copy result into X and Y
  for (int j=0; j<numColsRHS; j++) {
    w = dynamic_cast<NOX::LAPACK::Vector*>(&X[j]);
    for (int i=0; i<n; i++)
      (*w)(i) = RHS(i,j);
    for (int i=0; i<m; i++)
      Y(i,j) = RHS(n+i,j);
  }

  if (!res)
    return NOX::Abstract::Group::Failed;
  else
    return NOX::Abstract::Group::Ok;
  
}
