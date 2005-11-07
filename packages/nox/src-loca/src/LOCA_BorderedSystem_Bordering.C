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
  isContiguous(false)
{
}

LOCA::BorderedSystem::Bordering::~Bordering()
{
}

void
LOCA::BorderedSystem::Bordering::setIsContiguous(bool flag)
{
  isContiguous = flag;
}

void
LOCA::BorderedSystem::Bordering::setMatrixBlocks(
	 const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  grp = group;
  A = blockA;
  B = blockB;
  C = blockC;

  isZeroA = (A.get() == NULL);
  isZeroB = B->isDXZero();
  isZeroC = (C.get() == NULL);

  // ensure blocks B and C are not both zero
  if (isZeroB && isZeroC) 
    globalData->locaErrorCheck->throwError(
			    "LOCA::BorderedSystem::Bordering::setMatrixBlocks",
			    "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC) 
    globalData->locaErrorCheck->throwError(
			    "LOCA::BorderedSystem::Bordering::setMatrixBlocks",
			    "Blocks A and C cannot both be zero");
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

  isZeroF = (F == NULL);
  isZeroG = (G == NULL);

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

   if (isZeroA)
     status = solveAZero(params, B.get(), C.get(), F, G, X, Y);
   
   else if (isZeroB) {

     if (isContiguous) {
       Teuchos::RefCountPtr<NOX::Abstract::MultiVector> f = F->subView(indexF);
       Teuchos::RefCountPtr<NOX::Abstract::MultiVector> a = F->subView(indexA);
       Teuchos::RefCountPtr<NOX::Abstract::MultiVector> x = X.subView(indexF);

       status = solveBZero(params, a.get(), C.get(), f.get(), G, *x, Y);
     }
     else 
       status = solveBZero(params, A.get(), C.get(), F, G, X, Y);

   }
   
   else if (isZeroF)
     status = solveFZero(params, A.get(), B.get(), C.get(), G, X, Y);

   else {

     if (isContiguous) 
       status = solveContiguous(params, A.get(), B.get(), C.get(), 
				indexF, indexA, F, G, X, Y);

     else {
       int numColsRHS = numColsF + numColsA;
       Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS = 
	 F->clone(numColsRHS);
       Teuchos::RefCountPtr<NOX::Abstract::MultiVector> LHS = 
	 X.clone(numColsRHS);
       Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X1 = 
	 LHS->subView(indexF);
       RHS->setBlock(*F, indexF);
       RHS->setBlock(*A, indexA);
      
       status = solveContiguous(params, A.get(), B.get(), C.get(), 
				indexF, indexA, RHS.get(), G, *LHS, Y);
       X = *X1;
     }

   }

   return status;
}

// This function solves
//    | J A ||X|   |F|
//    | 0 C ||Y| = |G|
// via:  Y = C^-1 * G, X = J^-1 * (F - A*Y), where special cases of A,F,G=0
// are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveBZero(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::MultiVector* AA,
			    const NOX::Abstract::MultiVector::DenseMatrix* CC,
			    const NOX::Abstract::MultiVector* F,
			    const NOX::Abstract::MultiVector::DenseMatrix* G,
			    NOX::Abstract::MultiVector& X,
			    NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveBZero()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

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
	globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							       finalStatus,
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
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  else {
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS;

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
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}

// This function solves
//    | J   0 ||X|   |F|
//    | B^T C ||Y| = |G|
// via:  X = J^-1 * F, Y = C^-1 * (G - B^T*X), where special cases of B,F,G=0
// are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveAZero(
		       NOX::Parameter::List& params,
		       const LOCA::MultiContinuation::ConstraintInterface* BB,
		       const NOX::Abstract::MultiVector::DenseMatrix* CC,
		       const NOX::Abstract::MultiVector* F,
		       const NOX::Abstract::MultiVector::DenseMatrix* G,
		       NOX::Abstract::MultiVector& X,
		       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveAZero()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

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
    delete [] ipiv;
    if (info != 0) {
      status = NOX::Abstract::Group::Failed;
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							       finalStatus,
							       callingFunction);
    }
  }

  return finalStatus;
}

// This function solves
//    | J   A ||X|   |0|
//    | B^T C ||Y| = |G|
// via:  Xt = J^-1*A, Y = (C-B^T*Xt)^-1*G, X = -Xt*Y, where special cases of 
// C,G=0 are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveFZero(
		       NOX::Parameter::List& params,
		       const NOX::Abstract::MultiVector* AA,
		       const LOCA::MultiContinuation::ConstraintInterface* BB,
		       const NOX::Abstract::MultiVector::DenseMatrix* CC,
		       const NOX::Abstract::MultiVector::DenseMatrix* G,
		       NOX::Abstract::MultiVector& X,
		       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::Bordering::solveFZero()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Set X and Y to zero if G is zero
  if (isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return finalStatus;
  }

  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> Xt = 
    AA->clone(NOX::ShapeCopy);

  // compute Xt = J^-1 A
  status = grp->applyJacobianInverseMultiVector(params, *AA, *Xt);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  
  // compute t2 = -B^T*Xt
  NOX::Abstract::MultiVector::DenseMatrix t(BB->numConstraints(),
					    Xt->numVectors());
  BB->multiplyDX(-1.0, *Xt, t);
    
  // compute C - B^T*Xt
  if (!isZeroC)
    t += *CC;

  // compute Y = (C - B^T*Xt)^-1 * G
  Y.assign(*G);
  Teuchos::LAPACK<int,double> L;
  int *ipiv = new int[t.numRows()];
  int info;
  L.GESV(t.numRows(), Y.numCols(), t.values(), t.stride(), ipiv, 
	 Y.values(), Y.stride(), &info);
  delete [] ipiv;
  if (info != 0) {
    status = NOX::Abstract::Group::Failed;
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute X = -Xt*Y
  X.update(Teuchos::NO_TRANS, -1.0, *Xt, Y, -0.0);

  return finalStatus;
}

// This function assumes A, B, and F are nonzero.  C and/or G may be zero.
// It also assumes F and A are in a contiguous multivec, stored in F
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Bordering::solveContiguous(
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
    "LOCA::BorderedSystem::Bordering::solveContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // compute [X1 X2] = J^-1 [F A]
  status = grp->applyJacobianInverseMultiVector(params, *F, X);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X1 = X.subView(indexF);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X2 = X.subView(indexA);
  
  // compute t1 = -B^T*X1, for efficiency t1 is stored in Y
  BB->multiplyDX(-1.0, *X1, Y);
      
  // compute t2 = -B^T*X2
  NOX::Abstract::MultiVector::DenseMatrix t2(BB->numConstraints(),
					     X2->numVectors());
  BB->multiplyDX(-1.0, *X2, t2);

  // compute G - B^T*X1
  if (!isZeroG)
    Y += *G;
    
  // compute C - B^T*X2
  if (!isZeroC)
    t2 += *CC;

  // compute Y = (C - B^T*X2)^-1 * (G - B^T*X1)
  Teuchos::LAPACK<int,double> L;
  int *ipiv = new int[t2.numRows()];
  int info;
  L.GESV(t2.numRows(), Y.numCols(), t2.values(), t2.stride(), ipiv, 
	 Y.values(), Y.stride(), &info);
  delete [] ipiv;
  if (info != 0) {
    status = NOX::Abstract::Group::Failed;
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute X = X1 - X2*Y
  X1->update(Teuchos::NO_TRANS, -1.0, *X2, Y, 1.0);

  return finalStatus;
}
