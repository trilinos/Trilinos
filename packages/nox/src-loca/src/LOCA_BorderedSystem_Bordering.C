// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_BorderedSystem_Bordering.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_LAPACK.hpp"    // for LAPACK solve in applyInverse()

LOCA::BorderedSystem::Bordering::Bordering(NOX::Parameter::List& params): 
  grp(NULL),
  A(NULL),
  B(NULL),
  C(NULL),
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
  reset(params);
}

LOCA::BorderedSystem::Bordering::Bordering(
			       const LOCA::BorderedSystem::Bordering& source) :
  grp(source.grp),
  A(source.A),
  B(source.B),
  C(source.C),
  isZeroA(source.isZeroA),
  isZeroB(source.isZeroB),
  isZeroC(source.isZeroC),
  isZeroF(source.isZeroF),
  isZeroG(source.isZeroG),
  isContiguous(source.isContiguous),
  isZeroX(source.isZeroX),
  isZeroY(source.isZeroY),
  isZeroT1(source.isZeroT1),
  isZeroT2(source.isZeroT2)
{
}

LOCA::BorderedSystem::Bordering::~Bordering()
{
}

LOCA::BorderedSystem::Bordering&
LOCA::BorderedSystem::Bordering::operator=(
				const LOCA::BorderedSystem::Bordering& source)
{
  if (this != &source) {
    grp = source.grp;
    A = source.A;
    B = source.B;
    C = source.C;
    isZeroA = source.isZeroA;
    isZeroB = source.isZeroB;
    isZeroC = source.isZeroC;
    isZeroF = source.isZeroF;
    isZeroG = source.isZeroG;
    isContiguous = source.isContiguous;
    isZeroX = source.isZeroX;
    isZeroY = source.isZeroY;
    isZeroT1 = source.isZeroT1;
    isZeroT2 = source.isZeroT2;
  }
  
  return *this;
}

LOCA::BorderedSystem::Generic*
LOCA::BorderedSystem::Bordering::clone() const
{
  return new Bordering(*this);
}

LOCA::BorderedSystem::Generic& 
LOCA::BorderedSystem::Bordering::operator=(
				  const LOCA::BorderedSystem::Generic& source)
{
  return 
    operator=(dynamic_cast<const LOCA::BorderedSystem::Bordering&>(source));
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::reset(NOX::Parameter::List& params)
{
  return NOX::Abstract::Group::Ok;
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
    LOCA::ErrorCheck::throwError("LOCA::BorderedSystem::Bordering::setIsZero",
				 "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC) 
    LOCA::ErrorCheck::throwError("LOCA::BorderedSystem::Bordering::setIsZero",
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
     LOCA::ErrorCheck::throwError(
		     "LOCA::BorderedSystem::Bordering::setIsContiguous",
		     "Blocks F and A cannont be contiguous when one is zero");
}

void
LOCA::BorderedSystem::Bordering::setMatrixBlocks(
			const NOX::Abstract::Group* group,
			const NOX::Abstract::MultiVector* blockA,
			const NOX::Abstract::MultiVector* blockB,
			const NOX::Abstract::MultiVector::DenseMatrix* blockC)
{
  grp = group;
  A = blockA;
  B = blockB;
  C = blockC;
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
    X.multiply(1.0, *B, V);

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
    U.update(Teuchos::NO_TRANS, 1.0, *B, Y, 1.0);

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
   if (isContiguous) {
     for (int i=0; i<numColsA; i++)
       indexA[i] = numColsF + i;
   }
   else {
     for (int i=0; i<numColsA; i++)
       indexA[i] = i;
   }

   
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
       status = solveBZeroNoncontiguous(params, A, C, F, G, X, Y);

   }
   else {

     if (isContiguous  && !isZeroF && !isZeroA) 
       status = solveBNonZeroContiguous(params, A, B, C, indexF, indexA, 
					F, G, X, Y);

     else if (!isContiguous && !isZeroF && !isZeroA) {
       NOX::Abstract::MultiVector* RHS = F->clone(numColsRHS);
       NOX::Abstract::MultiVector* LHS = X.clone(numColsRHS);
       NOX::Abstract::MultiVector* X1 = LHS->subView(indexF);
       RHS->setBlock(*F, indexF);
       RHS->setBlock(*A, indexA);
      
       status = solveBNonZeroContiguous(params, A, B, C, indexF, indexA, RHS,
					G, *X1, Y);
       X = *X1;

       delete X1;
       delete RHS;
       delete LHS;
     }

     else if (isContiguous && (isZeroF || isZeroA)) {
       NOX::Abstract::MultiVector* f = F->subView(indexF);
       NOX::Abstract::MultiVector* a = F->subView(indexA);
       NOX::Abstract::MultiVector* x = X.subView(indexF);

       status = solveBNonZeroNoncontiguous(params, a, B, C, f, G, *x, Y);

       delete f;
       delete a; 
       delete x;
     }

     else 
       status = solveBNonZeroNoncontiguous(params, A, B, C, F, G, X, Y);
     
   }
   return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveBZeroNoncontiguous(
			      NOX::Parameter::List& params,
			      const NOX::Abstract::MultiVector* A,
			      const NOX::Abstract::MultiVector::DenseMatrix* C,
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
      RHS = A->clone(Y.numCols());
  }
    

  // Solve Y = C^-1 * G
  if (!isZeroY) {
    Y.assign(*G);
    M = new NOX::Abstract::MultiVector::DenseMatrix(*C);
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
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 1.0);
    else if (!isZeroA)
      RHS->update(Teuchos::NO_TRANS, -1.0, *A, Y, 0.0);

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
			      const NOX::Abstract::MultiVector* A,
			      const NOX::Abstract::MultiVector* B,
			      const NOX::Abstract::MultiVector::DenseMatrix* C,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveBZeroNoncontiguous()";
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
      X2 = A->clone(NOX::ShapeCopy);
      status = grp->applyJacobianInverseMultiVector(params, *A, *X2);
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }

  }

  if (!isZeroY) {

    // compute t1 = -B^T*X1, for efficiency t1 is stored in Y
    if (!isZeroT1) {
      X.multiply(-1.0, *B, Y);
    }

    // compute t2 = -B^T*X2
    if (!isZeroT2) {
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(B->numVectors(),
						       X2->numVectors());
      X2->multiply(-1.0, *B, *t2);
    }

    // compute G - B^T*X1
    if (!isZeroG && !isZeroT1)
      Y += *G;
    else if (!isZeroG)
      Y.assign(*G);  // don't use operator= (it destroys views)
    
    // compute C - B^T*X2
    if (!isZeroC && !isZeroT2)
      *t2 += *C;
    else if (!isZeroC)
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(*C);
      
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
			      const NOX::Abstract::MultiVector* A,
			      const NOX::Abstract::MultiVector* B,
			      const NOX::Abstract::MultiVector::DenseMatrix* C,
			      vector<int>& indexF,
			      vector<int>& indexA,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveBZeroContiguous()";
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
      X1->multiply(-1.0, *B, Y);
      
      // compute t2 = -B^T*X2
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(B->numVectors(),
						       X2->numVectors());
      X2->multiply(-1.0, *B, *t2);

    }

    // compute G - B^T*X1
    if (!isZeroG && !isZeroB)
      Y += *G;
    else if (!isZeroG)
      Y.assign(*G);  // don't use operator= (it destroys views)
    
    // compute C - B^T*X2
    if (!isZeroC && !isZeroB)
      *t2 += *C;
    else if (!isZeroC)
      t2 = new NOX::Abstract::MultiVector::DenseMatrix(*C);
      
    // Note that C and B cannot both be zero

    // compute Y = (C - B^T*X2)^-1 * (G - B^T*X1)
    Teuchos::LAPACK<int,double> L;
    int *ipiv = new int[t2->numRows()];
    int info;
    L.GESV(C->numRows(), Y.numCols(), t2->values(), t2->stride(), ipiv, 
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
