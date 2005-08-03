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

#include "LOCA_BorderedSystem_EpetraAugmented.H"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "NOX_Epetra_MultiVector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_EpetraNew_Group.H"
#include "LOCA_Epetra_AugmentedOp.H"

LOCA::BorderedSystem::EpetraAugmented::EpetraAugmented(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  A(),
  B(),
  C(),
  constraints(),
  numConstraints(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isContiguous(false)
{
}

LOCA::BorderedSystem::EpetraAugmented::~EpetraAugmented()
{
}

void
LOCA::BorderedSystem::EpetraAugmented::setIsContiguous(bool flag)
{
  isContiguous = flag;
}

void
LOCA::BorderedSystem::EpetraAugmented::setMatrixBlocks(
	 const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  string callingFunction = 
    "LOCA::BorderedSystem::EpetraAugmented::setMatrixBlocks";

  // Cast away const
  Teuchos::RefCountPtr<NOX::Abstract::Group> non_const_group = 
    Teuchos::rcp_const_cast<NOX::Abstract::Group>(group);

  // Cast group to an EpetraNew group
  grp = Teuchos::rcp_dynamic_cast<LOCA::EpetraNew::Group>(non_const_group);
  if (grp.get() == NULL)
    globalData->locaErrorCheck->throwError(
				    callingFunction,
				    "Group object must be an EpetraNew group");

  A = blockA;

  // Cast constraints to a ConstraintInterfaceMVDX
  constraints = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (constraints.get() == NULL)
    globalData->locaErrorCheck->throwError(
		 callingFunction,
		 "Constraints object must be of type ConstraintInterfaceMVDX");

  B = Teuchos::rcp(constraints->getDX(), false);
  C = blockC;

  // Determine which blocks are zero
  isZeroA = (A.get() == NULL);
  isZeroB = constraints->isDXZero();
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

  if (isZeroB)
    numConstraints = C->numRows();
  else
    numConstraints = B->numVectors();

  // We only use the augmented technique if A and B are nonzero 
  // (otherwise we use a block elimination).  If C is zero, we just create
  // a matrix of zeros and apply the standard algorithm
  if (isZeroC && !isZeroA && !isZeroB) {

    Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> tmpC = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							   B->numVectors(),
							   B->numVectors()));
    tmpC->putScalar(0.0);
    C = tmpC;
    isZeroC = false;
  }

}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::EpetraAugmented::apply(
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
LOCA::BorderedSystem::EpetraAugmented::applyTranspose(
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
LOCA::BorderedSystem::EpetraAugmented::applyInverse(
			      NOX::Parameter::List& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::EpetraAugmented::applyInverse()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  int numColsF;
  int numColsA;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F and G are zero, return zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  // ensure F and A are nonzero if contiguous
  if (isContiguous && (isZeroF || isZeroA)) 
     globalData->locaErrorCheck->throwError(
		    callingFunction,
		    "Blocks F and A cannont be contiguous when one is zero");

  if (!isZeroA)
     numColsA = A->numVectors();
   else
     numColsA = 0;

   if (!isZeroF && isContiguous)
     numColsF = F->numVectors() - numColsA;
   else if (!isZeroF)
     numColsF = F->numVectors();
   else
     numColsF = 0;

   // create subindexing vectors
   vector<int> indexF(numColsF);
   vector<int> indexA(numColsA);
   for (int i=0; i<numColsF; i++)
     indexF[i] = i;
   for (int i=0; i<numColsA; i++)
     indexA[i] = numColsF + i;
   
   Teuchos::RefCountPtr<NOX::Abstract::MultiVector> f;
   Teuchos::RefCountPtr<NOX::Abstract::MultiVector> a;
   Teuchos::RefCountPtr<NOX::Abstract::MultiVector> x;
   Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> cf;
   Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> ca;
   if (isContiguous) {
     f = Teuchos::rcp(F->subView(indexF));
     a = Teuchos::rcp(F->subView(indexA));
     x = Teuchos::rcp(X.subView(indexF));
     cf = f;
     ca = a;
   }
   else {
     cf = Teuchos::rcp(F, false);
     ca = A;
     x = Teuchos::rcp(&X, false);
   }

   if (isZeroA)
     finalStatus = solveAZero(params, constraints.get(), C.get(), F, G, X, Y);

   else if (isZeroB) 
     finalStatus = solveBZero(params, ca.get(), C.get(), cf.get(), G, *x, Y);

   else {
   
     // Get underlying Epetra vectors
     Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_a = 
       Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(ca);
     Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_b = 
       Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(B);
     Teuchos::RefCountPtr<NOX::Epetra::MultiVector> nox_epetra_x = 
       Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(x);
     Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_f = 
       Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(cf);
     Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_a = 
       Teuchos::rcp(&(nox_epetra_a->getEpetraMultiVector()), false);
     Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_b = 
       Teuchos::rcp(&(nox_epetra_b->getEpetraMultiVector()), false);
     Teuchos::RefCountPtr<Epetra_MultiVector> epetra_x = 
       Teuchos::rcp(&(nox_epetra_x->getEpetraMultiVector()), false);
     Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_f = 
       Teuchos::rcp(&(nox_epetra_f->getEpetraMultiVector()), false);
     
     // Get linear system
     NOX::EpetraNew::LinearSystem& linSys = grp->getLinearSystem();

     // Get Jacobian
     Teuchos::RefCountPtr<Epetra_Operator> jac =
       Teuchos::rcp(&linSys.getJacobianOperator(), false);

     // Set Jacobian
     linSys.setJacobianOperatorForSolve(*jac);

     // Create the preconditioner
     linSys.destroyPreconditioner();
     linSys.createPreconditioner(*((*epetra_x)(0)), params, false);

     // Get preconditioner
     Teuchos::RefCountPtr<Epetra_Operator> prec =
       Teuchos::rcp(&linSys.getGeneratedPrecOperator(), false);

     // Create augmented operators
     LOCA::Epetra::AugmentedOp extended_jac(globalData, jac, epetra_a, 
					    epetra_b, C);
     LOCA::Epetra::AugmentedOp extended_prec(globalData, prec, epetra_a, 
					     epetra_b, C);

     // Set augmented operator in linear system
     linSys.setJacobianOperatorForSolve(extended_jac);
     linSys.setPrecOperatorForSolve(extended_prec);

     // Create augmented Epetra vectors for x, f
     Teuchos::RefCountPtr<Epetra_MultiVector> epetra_augmented_x = 
       extended_jac.buildEpetraAugmentedMultiVec(*epetra_x, &Y, false);
     Teuchos::RefCountPtr<Epetra_MultiVector> epetra_augmented_f = 
       extended_jac.buildEpetraAugmentedMultiVec(*epetra_f, G, true);

     // Create augmented NOX::Epetra::MultiVectors as views
     NOX::Epetra::MultiVector nox_epetra_augmented_x(*epetra_augmented_x,
						     NOX::DeepCopy, true);
     NOX::Epetra::MultiVector nox_epetra_augmented_f(*epetra_augmented_f,
						     NOX::DeepCopy, true);

     // Solve for each RHS
     int m = nox_epetra_augmented_x.numVectors();
     for (int i=0; i<m; i++) {
       extended_jac.init(*((*epetra_f)(i)));
       extended_prec.init(*((*epetra_f)(i)));
       bool stat = 
	 linSys.applyJacobianInverse(
		params, 
		dynamic_cast<NOX::Epetra::Vector&>(nox_epetra_augmented_f[i]),
		dynamic_cast<NOX::Epetra::Vector&>(nox_epetra_augmented_x[i]));
       if (stat == true)
	 status = NOX::Abstract::Group::Ok;
       else
	 status = NOX::Abstract::Group::NotConverged;
       finalStatus = 
	 LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						      callingFunction);
     }

     // Set results
     extended_jac.init(*epetra_x);
     extended_jac.setEpetraAugmentedMultiVec(*epetra_x, Y,
					     *epetra_augmented_x);

     // Set original Jacobian in linear system
     linSys.setJacobianOperatorForSolve(*jac);
     linSys.destroyPreconditioner();
     //linSys.setPrecOperatorForSolve(*prec);
   }

   return finalStatus;

}

// This function solves
//    | J   0 ||X|   |F|
//    | B^T C ||Y| = |G|
// via:  X = J^-1 * F, Y = C^-1 * (G - B^T*X), where special cases of B,F,G=0
// are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::EpetraAugmented::solveAZero(
		       NOX::Parameter::List& params,
		       const LOCA::MultiContinuation::ConstraintInterface* BB,
		       const NOX::Abstract::MultiVector::DenseMatrix* CC,
		       const NOX::Abstract::MultiVector* F,
		       const NOX::Abstract::MultiVector::DenseMatrix* G,
		       NOX::Abstract::MultiVector& X,
		       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::EpetraAugmented::solveAZero()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);


  // Determine if X or Y is zero
  bool isZeroX = isZeroF;
  bool isZeroY = isZeroG && (isZeroB || isZeroX);

  // First compute X
  if (isZeroX)
    X.init(0.0);
  else {
    // Solve X = J^-1 F, note F must be nonzero
    status = grp->applyJacobianInverseMultiVector(params, *F, X);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Now compute Y
  if (isZeroY)
    Y.putScalar(0.0);
  else {
    // Compute G - B^T*X and store in Y
    if (isZeroG) 
      BB->multiplyDX(-1.0, X, Y);
    else {
      Y.assign(*G);
      if (!isZeroB && !isZeroX) {
	NOX::Abstract::MultiVector::DenseMatrix T(Y.numRows(),Y.numCols());
	BB->multiplyDX(1.0, X, T);
	Y -= T;
      }
    }

    // Overwrite Y with Y = C^-1 * (G - B^T*X)
    NOX::Abstract::MultiVector::DenseMatrix M(*CC);
    int *ipiv = new int[M.numRows()];
    Teuchos::LAPACK<int,double> L;
    int info;
    L.GESV(M.numRows(), Y.numCols(), M.values(), M.stride(), ipiv, 
	   Y.values(), Y.stride(), &info);
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }
    delete [] ipiv;
  }

  return finalStatus;
}

// This function solves
//    | J A ||X|   |F|
//    | 0 C ||Y| = |G|
// via:  Y = C^-1 * G, X = J^-1 * (F - A*Y), where special cases of A,F,G=0
// are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::EpetraAugmented::solveBZero(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::MultiVector* AA,
			    const NOX::Abstract::MultiVector::DenseMatrix* CC,
			    const NOX::Abstract::MultiVector* F,
			    const NOX::Abstract::MultiVector::DenseMatrix* G,
			    NOX::Abstract::MultiVector& X,
			    NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::EpetraAugmented::solveBZero()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // Determine if X or Y is zero
  bool isZeroY = isZeroG;
  bool isZeroX = isZeroF && (isZeroA || isZeroY);

  // First compute Y
  if (isZeroY)
    Y.putScalar(0.0);
  else {
    // Solve Y = C^-1 * G
    NOX::Abstract::MultiVector::DenseMatrix M(*CC);
    int *ipiv = new int[M.numRows()];
    Teuchos::LAPACK<int,double> L;
    int info;
    
    Y.assign(*G);
    L.GESV(M.numRows(), Y.numCols(), M.values(), M.stride(), ipiv, 
	   Y.values(), Y.stride(), &info);
    delete [] ipiv;
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }
  }

  // Now compute X
  if (isZeroX)
    X.init(0.0);
  else if (isZeroA || isZeroY) {
    // Solve X = J^-1 F, note F must be nonzero
    status = grp->applyJacobianInverseMultiVector(params, *F, X);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  else {
    NOX::Abstract::MultiVector *RHS;

    if (isZeroF) {
      RHS = AA->clone(Y.numCols());
      RHS->update(Teuchos::NO_TRANS, -1.0, *AA, Y, 0.0);
    }
    else {
      RHS = F->clone(NOX::DeepCopy);
      RHS->update(Teuchos::NO_TRANS, -1.0, *AA, Y, 1.0);
    }
    // Solve X = J^-1 (F-A*Y)
    status = grp->applyJacobianInverseMultiVector(params, *RHS, X);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
    delete RHS;
  }

  return finalStatus;
}
