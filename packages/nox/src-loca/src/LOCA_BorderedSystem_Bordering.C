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
  isZeroC(true)
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
  isZeroC(source.isZeroC)
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
LOCA::BorderedSystem::Bordering::setIsZero(bool flagA, bool flagB, bool flagC)
{
  isZeroA = flagA;
  isZeroB = flagB;
  isZeroC = flagC;

  // ensure blocks B and C are not both zero
  if (isZeroB && isZeroC) 
    LOCA::ErrorCheck::throwError("LOCA::BorderedSystem::Bordering::setIsZero",
				 "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC) 
    LOCA::ErrorCheck::throwError("LOCA::BorderedSystem::Bordering::setIsZero",
				 "Blocks A and C cannot both be zero");
}

void
LOCA::BorderedSystem::Bordering::setBlocks(
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
			    bool isZeroF,
			    bool isZeroG,
			    bool contiguousRHS,
			    const NOX::Abstract::MultiVector* F,
			    const NOX::Abstract::MultiVector::DenseMatrix* G,
			    NOX::Abstract::MultiVector& X,
		            NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::applyInverse()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  NOX::Abstract::MultiVector* LHS = NULL;
  NOX::Abstract::MultiVector* RHS = NULL;
  const NOX::Abstract::MultiVector* cRHS = NULL;
  NOX::Abstract::MultiVector* X1 = NULL;
  NOX::Abstract::MultiVector* X2 = NULL;
  NOX::Abstract::MultiVector::DenseMatrix* t2 = NULL;
  const NOX::Abstract::MultiVector::DenseMatrix* t3 = NULL;

  bool isZeroX = isZeroF && isZeroA;
  bool isZeroY = (isZeroG && isZeroB) || (isZeroG && isZeroF);
  bool isZeroT1 = isZeroB || isZeroF;
  bool isZeroT2 = isZeroB || isZeroA;

  // set X and Y if they are zero
  if (isZeroX)
    X.init(0.0);

  if (isZeroY)
    Y.putScalar(0.0);

  if (!isZeroX) {

    // form concatenated right-hand and left-hand sides
    if (!isZeroF && !isZeroA) {

      int numColsF;
      int numColsA;
      int numColsRHS;

      // get number of columns for F and A
      if (!contiguousRHS) {
	numColsF = F->numVectors();
	numColsA = A->numVectors();
	numColsRHS = numColsF + numColsA;
      }
      else {
	numColsA = A->numVectors();
	numColsF = F->numVectors() - numColsA;
	numColsRHS = numColsF + numColsA;
      }

      // create subindexing vectors
      vector<int> indexF(numColsF);
      vector<int> indexA(numColsA);
      for (int i=0; i<numColsF; i++)
	indexF[i] = i;
      for (int i=0; i<numColsA; i++)
	indexA[i] = numColsF + i;

      // concatenate F and A into one multivec if they aren't contiguous
      if (!contiguousRHS) {
	RHS = F->clone(numColsRHS);
	RHS->setBlock(*F, indexF);
	RHS->setBlock(*A, indexA);
	cRHS = RHS;

	LHS = X.clone(numColsRHS);
      }
      else  {
	cRHS = F;
	LHS = &X;
      }

      // compute J^-1 [F A]
      status = grp->applyJacobianInverseMultiVector(params, *cRHS, *LHS);
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);

      // Get views for J^-1 F and J^-1 A
      X1 = LHS->subView(indexF);
      X2 = LHS->subView(indexA);
    }

    else {

      X1 = &X;

      // compute X1 = J^-1*F
      if (!isZeroF) {
	status = grp->applyJacobianInverseMultiVector(params, *F, *X1);
	finalStatus = 
	  LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						       callingFunction);
      }

      // compute X2 = J^-1*A
      if (!isZeroA) {
	X2 = A->clone(NOX::ShapeCopy);
	status = grp->applyJacobianInverseMultiVector(params, *A, *X2);
	finalStatus = 
	  LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						       callingFunction);
      }

    }

  }

  if (!isZeroY) {

    // compute t1 = -B^T*X1, for efficiency t1 is stored in Y
    if (!isZeroT1) {
      X1->multiply(-1.0, *B, Y);
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
    if (!isZeroC && !isZeroT2) {
      *t2 += *C;
      t3 = t2;
    }
    else if (!isZeroC)
      t3 = C;

    // Note that C and T2 cannot both be zero

    // compute Y = (C - B^T*X2)^-1 * (G - B^T*X1)
    Teuchos::LAPACK<int,double> L;
    int info;
    int *ipiv = new int[C->numRows()];
    L.GESV(C->numRows(), Y.numCols(), t3->values(), t3->stride(), ipiv, 
	   Y.values(), Y.stride(), &info);
    delete [] ipiv;
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus = 
	LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						     callingFunction);
    }

  }

  if (!isZeroX) {
    
    // compute X = X1 - X2*Y

    // if A is zero or Y is zero, then nothing needs to be done
    if (!isZeroA && !isZeroY) {
      if (isZeroF)
	X1->update(Teuchos::NO_TRANS, -1.0, *X2, Y, 0.0);

      else
	X1->update(Teuchos::NO_TRANS, -1.0, *X2, Y, 1.0);
    }

    // copy X1 into X in the one case when they are not equal
    if (!isZeroF && !isZeroA && !contiguousRHS) 
      X = *X1;
    
  }

  // clean up temporaries
  if (!isZeroF && !isZeroA && !contiguousRHS) {
    delete LHS;
    delete RHS;
  }
  if (!isZeroF && !isZeroA && contiguousRHS)
    delete X1;
  if (!isZeroA)
    delete X2;
  if (!isZeroT2)
    delete t2;

  return finalStatus;
}
