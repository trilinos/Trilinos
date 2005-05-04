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

#include "LOCA_BorderedSystem_LAPACKDirectSolve.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_LAPACK.H"

LOCA::BorderedSystem::LAPACKDirectSolve::LAPACKDirectSolve(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  A(),
  B(),
  C(),
  augmentedJ(),
  pivots(),
  n(0),
  m(0),
  N(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isZeroF(true),
  isZeroG(true),
  isContiguous(false)
{
}

LOCA::BorderedSystem::LAPACKDirectSolve::~LAPACKDirectSolve()
{
}

void
LOCA::BorderedSystem::LAPACKDirectSolve::setIsZero(bool flagA, bool flagB, 
						   bool flagC,
						   bool flagF, bool flagG)
{
  isZeroA = flagA;
  isZeroB = flagB;
  isZeroC = flagC;
  isZeroF = flagF;
  isZeroG = flagG;

  // ensure blocks B and C are not both zero
  if (isZeroB && isZeroC) 
    globalData->locaErrorCheck->throwError(
			"LOCA::BorderedSystem::LAPACKDirectSolve::setIsZero",
			"Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC) 
    globalData->locaErrorCheck->throwError(
			"LOCA::BorderedSystem::LAPACKDirectSolve::setIsZero",
			"Blocks A and C cannot both be zero");

}

void
LOCA::BorderedSystem::LAPACKDirectSolve::setIsContiguous(bool flag)
{
  isContiguous = flag;

  // ensure F and A are nonzero if contiguous
  if (isContiguous && (isZeroF || isZeroA)) 
     globalData->locaErrorCheck->throwError(
		  "LOCA::BorderedSystem::LAPACKDirectSolve::setIsContiguous",
		  "Blocks F and A cannont be contiguous when one is zero");
}

void
LOCA::BorderedSystem::LAPACKDirectSolve::setMatrixBlocks(
	 const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockBC)
{
  string callingFunction = 
    "LOCA::BorderedSystem::LAPACKDirectSolve::setMatrixBlocks()";
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
  B = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockBC);
  if (B.get() == NULL)
    globalData->locaErrorCheck->throwError(
	     callingFunction,
	     string("Constraint argument is not of type\n") + 
	     string("LOCA::MultiContinuation::ConstraintInterfaceMVDX!\n") +
	     string("The LAPACK Direct Solve bordered solver method can\n") +
	     string("only be used with constraints that support obtaining\n") +
	     string("the constraint derivative as a multivector."));
    
  C = B->getConstraintDerivativesP();

  // Get the Jacobian matrix and size
  const NOX::LAPACK::Matrix& J = grp->getJacobianMatrix();
  n = J.numRows();

  // Get the number of additional rows/columns
  m = A->numVectors();

  // Form a new (n+m) x (n+m) matrix if this is a new size
  if (n+m != N) {
    N = n+m;
    augmentedJ = Teuchos::rcp(new NOX::LAPACK::Matrix(N,N));
    pivots.resize(N);
  }

  // Copy Jacobian
  for (int j=0; j<n; j++)
    for (int i=0; i<n; i++)
      (*augmentedJ)(i,j) = J(i,j);

  // Copy A
  if (isZeroA) {
    for (int j=0; j<m; j++) 
      for (int i=0; i<n; i++)
	(*augmentedJ)(i,j+n) = 0.0;
  }
  else {
    for (int j=0; j<m; j++) {
      v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*A)[j]);
      for (int i=0; i<n; i++)
	(*augmentedJ)(i,j+n) = (*v)(i);
    }
  }

  // Copy B
  if (isZeroB) {
    for (int i=0; i<m; i++) 
      for (int j=0; j<n; j++)
	(*augmentedJ)(i+n,j) = 0.0;
  }
  else {
    BV = B->getConstraintDerivativesX();
    for (int i=0; i<m; i++) {
      v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*BV)[i]);
      for (int j=0; j<n; j++)
	(*augmentedJ)(i+n,j) = (*v)(j);
    }
  }

  // Copy C
  if (isZeroC) {
    for (int j=0; j<m; j++)
      for (int i=0; i<m; i++)
	(*augmentedJ)(i+n,j+n) = 0.0;
  }
  else {
    for (int j=0; j<m; j++)
      for (int i=0; i<m; i++)
	(*augmentedJ)(i+n,j+n) = (*C)(i,j);
  }

  // Factor augmented Jacobian
  int info;
  DGETRF_F77(&N, &N, &(*augmentedJ)(0,0), &N, &pivots[0], &info);
  if (info != 0)
    globalData->locaErrorCheck->throwError(
		 "LOCA::borderedSystem::LAPACKDirectSolve::setMatrixBlocks()",
		 "Factorization of augmented Jacobian matrix failed!");
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::LAPACKDirectSolve::apply(
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
    B->applyConstraintDerivativesX(1.0, X, V);

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
LOCA::BorderedSystem::LAPACKDirectSolve::applyTranspose(
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
    B->applyConstraintDerivativesX(Teuchos::NO_TRANS, 1.0, Y, 1.0, U);

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
LOCA::BorderedSystem::LAPACKDirectSolve::applyInverse(
			      NOX::Parameter::List& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  // If F & G are zero, the solution is zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  int numColsRHS;
  const NOX::LAPACK::Vector *v;
  NOX::LAPACK::Vector *w;

  if (!isZeroF && isContiguous)
    numColsRHS = F->numVectors() - A->numVectors();
  else if (!isZeroF)
    numColsRHS = F->numVectors();
  else 
    numColsRHS = G->numCols();

  const NOX::Abstract::MultiVector *FF;
  if (!isZeroF && isContiguous) {
    // create subindexing vectors
    vector<int> indexF(numColsRHS);
    for (int i=0; i<numColsRHS; i++)
      indexF[i] = i;
    
    // Get actual F
    FF = F->subView(indexF);
  }
  else if (!isZeroF)
    FF = F;
  else
    FF = NULL;
    
  // Concatenate F & G into a single matrix
  NOX::LAPACK::Matrix RHS(N,numColsRHS);
  if (isZeroF) {
    for (int j=0; j<numColsRHS; j++)
      for (int i=0; i<n; i++)
	RHS(i,j) = 0.0;
  }
  else {
    for (int j=0; j<numColsRHS; j++) {
      v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*FF)[j]);
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
  int info;
  DGETRS_F77("N", &N, &numColsRHS, &(*augmentedJ)(0,0), &N, &pivots[0],
	     &RHS(0,0), &N, &info);

  // Copy result into X and Y
  for (int j=0; j<numColsRHS; j++) {
    w = dynamic_cast<NOX::LAPACK::Vector*>(&X[j]);
    for (int i=0; i<n; i++)
      (*w)(i) = RHS(i,j);
    for (int i=0; i<m; i++)
      Y(i,j) = RHS(n+i,j);
  }

  if (!isZeroF && isContiguous)
    delete FF;

  if (info != 0)
    return NOX::Abstract::Group::Failed;
  else
    return NOX::Abstract::Group::Ok;
  
}
