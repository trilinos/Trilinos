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

#include "LOCA_BorderedSystem_Bordering.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "Teuchos_LAPACK.hpp"    // for LAPACK solve in applyInverse()

LOCA::BorderedSystem::Bordering::Bordering(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  A(),
  B(),
  C(),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isZeroF(true),
  isZeroG(true),
  isContiguous(false),
  isZeroX(true),
  isZeroY(true),
  isZeroT1(true),
  isZeroT2(true)
{
}

LOCA::BorderedSystem::Bordering::~Bordering()
{
}

void
LOCA::BorderedSystem::Bordering::setIsZero(bool flagA, bool flagB, bool flagC,
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
				 "LOCA::BorderedSystem::Bordering::setIsZero",
				 "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC) 
    globalData->locaErrorCheck->throwError(
				 "LOCA::BorderedSystem::Bordering::setIsZero",
				 "Blocks A and C cannot both be zero");

  isZeroX = isZeroF && isZeroA;
  isZeroY = (isZeroG && isZeroB) || (isZeroG && isZeroF);
  isZeroT1 = isZeroB || isZeroF;
  isZeroT2 = isZeroB || isZeroA;
}

void
LOCA::BorderedSystem::Bordering::setIsContiguous(bool flag)
{
  isContiguous = flag;

  // ensure F and A are nonzero if contiguous
  if (isContiguous && (isZeroF || isZeroA)) 
     globalData->locaErrorCheck->throwError(
		     "LOCA::BorderedSystem::Bordering::setIsContiguous",
		     "Blocks F and A cannont be contiguous when one is zero");
}

void
LOCA::BorderedSystem::Bordering::setMatrixBlocks(
	 const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockBC)
{
  grp = group;
  A = blockA;
  B = blockBC;
  C = &(B->getConstraintDerivativesP());
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::apply(
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
LOCA::BorderedSystem::Bordering::applyTranspose(
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
LOCA::BorderedSystem::Bordering::applyInverse(
			      NOX::Parameter::List& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::applyInverse()";
  NOX::Abstract::Group::ReturnType status;

  int numColsF;
  int numColsA;
  int numColsRHS;

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

   numColsRHS = numColsF + numColsA;

    // create subindexing vectors
   vector<int> indexF(numColsF);
   vector<int> indexA(numColsA);
   for (int i=0; i<numColsF; i++)
     indexF[i] = i;
   for (int i=0; i<numColsA; i++)
     indexA[i] = numColsF + i;

   
   if (isZeroB) {

     if (isContiguous) {
       NOX::Abstract::MultiVector* f = F->subView(indexF);
       NOX::Abstract::MultiVector* a = F->subView(indexA);
       NOX::Abstract::MultiVector* x = X.subView(indexF);

       status = solveBZeroNoncontiguous(params, a, C, f, G, *x, Y);

       delete f;
       delete a; 
       delete x;
     }
     else 
       status = solveBZeroNoncontiguous(params, A.get(), C, 
					F, G, X, Y);

   }
   else {

     if (isContiguous  && !isZeroF && !isZeroA) 
       status = solveBNonZeroContiguous(params, A.get(), B.get(), C, 
					indexF, indexA, 
					F, G, X, Y);

     else if (!isContiguous && !isZeroF && !isZeroA) {
       NOX::Abstract::MultiVector* RHS = F->clone(numColsRHS);
       NOX::Abstract::MultiVector* LHS = X.clone(numColsRHS);
       NOX::Abstract::MultiVector* X1 = LHS->subView(indexF);
       RHS->setBlock(*F, indexF);
       RHS->setBlock(*A, indexA);
      
       status = solveBNonZeroContiguous(params, A.get(), B.get(), C, 
					indexF, indexA, RHS,
					G, *LHS, Y);
       X = *X1;

       delete X1;
       delete RHS;
       delete LHS;
     }

     else if (isContiguous && (isZeroF || isZeroA)) {
       NOX::Abstract::MultiVector* f = F->subView(indexF);
       NOX::Abstract::MultiVector* a = F->subView(indexA);
       NOX::Abstract::MultiVector* x = X.subView(indexF);

       status = solveBNonZeroNoncontiguous(params, a, B.get(), C, f, 
					   G, *x, Y);

       delete f;
       delete a; 
       delete x;
     }

     else 
       status = solveBNonZeroNoncontiguous(params, A.get(), B.get(), 
					   C, F, G, X, Y);
     
   }
   return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveBZeroNoncontiguous(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::MultiVector* AA,
			    const NOX::Abstract::MultiVector::DenseMatrix* CC,
			    const NOX::Abstract::MultiVector* F,
			    const NOX::Abstract::MultiVector::DenseMatrix* G,
			    NOX::Abstract::MultiVector& X,
			    NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveBZeroNoncontiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  NOX::Abstract::MultiVector *RHS;
  NOX::Abstract::MultiVector::DenseMatrix *M;
  Teuchos::LAPACK<int,double> L;
  int *ipiv;
  int info;

  // set X and Y if they are zero
  if (isZeroX)
    X.init(0.0);

  if (isZeroY)
    Y.putScalar(0.0);

  if (!isZeroX) {
    if (!isZeroF)
      RHS = F->clone(NOX::DeepCopy);
    else
      RHS = AA->clone(Y.numCols());
  }
    

  // Solve Y = C^-1 * G
  if (!isZeroY) {
    Y.assign(*G);
    M = new NOX::Abstract::MultiVector::DenseMatrix(*CC);
    ipiv = new int[M->numRows()];
    
    L.GESV(M->numRows(), Y.numCols(), M->values(), M->stride(), ipiv, 
	   Y.values(), Y.stride(), &info);
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }
    delete [] ipiv;
    delete M;

    // compute F - A*Y
    if (!isZeroF && !isZeroA)
      RHS->update(Teuchos::NO_TRANS, -1.0, *AA, Y, 1.0);
    else if (!isZeroA)
      RHS->update(Teuchos::NO_TRANS, -1.0, *AA, Y, 0.0);

  }

  // Solve X = J^-1 (F-A*Y)
  if (!isZeroX) {
    status = grp->applyJacobianInverseMultiVector(params, *RHS, X);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
    delete RHS;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveBNonZeroNoncontiguous(
		      NOX::Parameter::List& params,
		      const NOX::Abstract::MultiVector* AA,
		      const LOCA::MultiContinuation::ConstraintInterface* BB,
		      const NOX::Abstract::MultiVector::DenseMatrix* CC,
		      const NOX::Abstract::MultiVector* F,
		      const NOX::Abstract::MultiVector::DenseMatrix* G,
		      NOX::Abstract::MultiVector& X,
		      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveBNonZeroNoncontiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // set X and Y if they are zero
  if (isZeroX)
    X.init(0.0);

  if (isZeroY)
    Y.putScalar(0.0);

  NOX::Abstract::MultiVector *X2;
  NOX::Abstract::MultiVector::DenseMatrix *t2;

  if (!isZeroX) {

    // compute X1 = J^-1*F, X1 is stored in X
    if (!isZeroF) {
      status = grp->applyJacobianInverseMultiVector(params, *F, X);
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }

    // compute X2 = J^-1*A
    if (!isZeroA && !isZeroY) {
      X2 = AA->clone(NOX::ShapeCopy);
      status = grp->applyJacobianInverseMultiVector(params, *AA, *X2);
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }

  }

  if (!isZeroY) {

    // compute t1 = -B^T*X1, for efficiency t1 is stored in Y
    if (!isZeroT1) {
      BB->applyConstraintDerivativesX(-1.0, X, Y);
    }

    // compute t2 = -B^T*X2
    if (!isZeroT2) {
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(BB->numConstraints(),
						       X2->numVectors());
      BB->applyConstraintDerivativesX(-1.0, *X2, *t2);
    }

    // compute G - B^T*X1
    if (!isZeroG && !isZeroT1)
      Y += *G;
    else if (!isZeroG)
      Y.assign(*G);  // don't use operator= (it destroys views)
    
    // compute C - B^T*X2
    if (!isZeroC && !isZeroT2)
      *t2 += *CC;
    else if (!isZeroC)
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(*CC);
      
    // Note that C and T2 cannot both be zero

    // compute Y = (C - B^T*X2)^-1 * (G - B^T*X1)
    Teuchos::LAPACK<int,double> L;
    int *ipiv = new int[t2->numRows()];
    int info;
    L.GESV(t2->numRows(), Y.numCols(), t2->values(), t2->stride(), ipiv, 
	   Y.values(), Y.stride(), &info);
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }
    delete [] ipiv;
    delete t2;

  }

  if (!isZeroX) {

    // compute X = X1 - X2*Y
      
    // if A is zero or Y is zero, then nothing needs to be done
    if (!isZeroA && !isZeroY) {
      if (isZeroF)
	X.update(Teuchos::NO_TRANS, -1.0, *X2, Y, 0.0);
      
      else
	X.update(Teuchos::NO_TRANS, -1.0, *X2, Y, 1.0);

      delete X2;
    }

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveBNonZeroContiguous(
		       NOX::Parameter::List& params,
		       const NOX::Abstract::MultiVector* AA,
		       const LOCA::MultiContinuation::ConstraintInterface* BB,
		       const NOX::Abstract::MultiVector::DenseMatrix* CC,
		       vector<int>& indexF,
		       vector<int>& indexA,
		       const NOX::Abstract::MultiVector* F,
		       const NOX::Abstract::MultiVector::DenseMatrix* G,
		       NOX::Abstract::MultiVector& X,
		       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveBNonZeroContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  NOX::Abstract::MultiVector *X1;
  NOX::Abstract::MultiVector *X2; 
  NOX::Abstract::MultiVector::DenseMatrix *t2;
  
  // set X and Y if they are zero
  if (isZeroX)
    X.init(0.0);

  if (isZeroY)
    Y.putScalar(0.0);

  if (!isZeroX) {

    // compute J^-1 [F A]
    status = grp->applyJacobianInverseMultiVector(params, *F, X);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    X1 = X.subView(indexF);
    X2 = X.subView(indexA);
  }

  if (!isZeroY) {

    
    if (!isZeroB) {

      // compute t1 = -B^T*X1, for efficiency t1 is stored in Y
      BB->applyConstraintDerivativesX(-1.0, *X1, Y);
      
      // compute t2 = -B^T*X2
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(BB->numConstraints(),
						       X2->numVectors());
      BB->applyConstraintDerivativesX(-1.0, *X2, *t2);

    }

    // compute G - B^T*X1
    if (!isZeroG && !isZeroB)
      Y += *G;
    else if (!isZeroG)
      Y.assign(*G);  // don't use operator= (it destroys views)
    
    // compute C - B^T*X2
    if (!isZeroC && !isZeroB)
      *t2 += *CC;
    else if (!isZeroC)
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(*CC);
      
    // Note that C and B cannot both be zero

    // compute Y = (C - B^T*X2)^-1 * (G - B^T*X1)
    Teuchos::LAPACK<int,double> L;
    int *ipiv = new int[t2->numRows()];
    int info;
    L.GESV(CC->numRows(), Y.numCols(), t2->values(), t2->stride(), ipiv, 
	   Y.values(), Y.stride(), &info);
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }
    delete t2;

  }

  // compute X = X1 - X2*Y
  if (!isZeroX) {
    X1->update(Teuchos::NO_TRANS, -1.0, *X2, Y, 1.0);
    delete X1;
    delete X2;
  }

  return finalStatus;
}
