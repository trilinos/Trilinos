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

#include "LOCA_BorderedSolver_Bordering.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "Teuchos_LAPACK.hpp"    // for LAPACK solve in applyInverse()
#include "LOCA_BorderedSolver_LowerTriangularBlockElimination.H"
#include "LOCA_BorderedSolver_UpperTriangularBlockElimination.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"

LOCA::BorderedSolver::Bordering::Bordering(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  ts_grp(),
  A(),
  B(),
  C(),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isZeroF(true),
  isZeroG(true)
{
}

LOCA::BorderedSolver::Bordering::~Bordering()
{
}

void
LOCA::BorderedSolver::Bordering::setMatrixBlocks(
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
			    "LOCA::BorderedSolver::Bordering::setMatrixBlocks",
			    "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC) 
    globalData->locaErrorCheck->throwError(
			    "LOCA::BorderedSolver::Bordering::setMatrixBlocks",
			    "Blocks A and C cannot both be zero");

  // cast group to a transpose-solve group for applyInverseTranspose()
  ts_grp = 
    Teuchos::rcp_dynamic_cast<const LOCA::Abstract::TransposeSolveGroup>(grp);
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Bordering::initForSolve()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Bordering::initForTransposeSolve()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Bordering::apply(
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
LOCA::BorderedSolver::Bordering::applyTranspose(
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
LOCA::BorderedSolver::Bordering::applyInverse(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSolver::Bordering::applyInverse()";
  NOX::Abstract::Group::ReturnType status;

  isZeroF = (F == NULL);
  isZeroG = (G == NULL);

   if (isZeroA) {
     LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
     status = ltbe.solve(params, *grp, *B, *C, F, G, X, Y);
   }
   
   else if (isZeroB) {
     LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
     status = utbe.solve(params, *grp, A.get(), *C, F, G, X, Y);

   }
   
   else if (isZeroF)
     status = solveFZero(params, A.get(), B.get(), C.get(), G, X, Y);

   else {

     int numColsA = A->numVectors();
     int numColsF = F->numVectors();

     // create indexing vectors
     vector<int> indexF(numColsF);
     vector<int> indexA(numColsA);
     for (int i=0; i<numColsF; i++)
       indexF[i] = i;
     for (int i=0; i<numColsA; i++)
       indexA[i] = numColsF + i;
     int numColsRHS = numColsF + numColsA;

     // copy F & A into 1 multivector
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS = 
       F->clone(numColsRHS);
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> LHS = 
       X.clone(numColsRHS);
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X1 = 
       LHS->subView(indexF);
     RHS->setBlock(*F, indexF);
     RHS->setBlock(*A, indexA);
      
     // solve
     status = solveContiguous(params, A.get(), B.get(), C.get(), 
			      indexF, indexA, RHS.get(), G, *LHS, Y);
     X = *X1;

   }

   return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Bordering::applyInverseTranspose(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSolver::Bordering::applyInverseTranspose()";
  NOX::Abstract::Group::ReturnType status;

  isZeroF = (F == NULL);
  isZeroG = (G == NULL);

  // The group must implement the transpose solve interface 
  if (ts_grp == Teuchos::null)
    globalData->locaErrorCheck->throwError(
	 callingFunction,
	 string("Group must implement the LOCA::Abstract::TransposeSolveGroup")
	 + string(" interface in order to solve the transpose of the bordered")
	 + string(" system via bordering."));

  // For the transpose solve, B must be a multi-vec constraint if it is nonzero
  Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> B_mvdx;
  const NOX::Abstract::MultiVector* BB = NULL;

  if (!isZeroB) {
    B_mvdx = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(B);
    if (B == Teuchos::null)
      globalData->locaErrorCheck->throwError(
		 callingFunction,
		 "Constraints object must be of type ConstraintInterfaceMVDX");
    BB = B_mvdx->getDX();
  }

   if (isZeroA) {
     LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
     status = utbe.solveTranspose(params, *ts_grp, BB, *C, F, G, X, Y);
   }
   
   else if (isZeroB) {
     LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
     status = ltbe.solveTranspose(params, *ts_grp, *A, *C, F, G, X, Y);

   }
   
   else if (isZeroF)
     status = solveFZeroTrans(params, A.get(), BB, C.get(), G, X, Y);

   else {

     int numColsB = BB->numVectors();
     int numColsF = F->numVectors();

     // create indexing vectors
     vector<int> indexF(numColsF);
     vector<int> indexB(numColsB);
     for (int i=0; i<numColsF; i++)
       indexF[i] = i;
     for (int i=0; i<numColsB; i++)
       indexB[i] = numColsF + i;
     int numColsRHS = numColsF + numColsB;

     // copy F & A into 1 multivector
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS = 
       F->clone(numColsRHS);
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> LHS = 
       X.clone(numColsRHS);
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X1 = 
       LHS->subView(indexF);
     RHS->setBlock(*F, indexF);
     RHS->setBlock(*BB, indexB);
      
     // solve
     status = solveContiguousTrans(params, A.get(), BB, C.get(), 
				   indexF, indexB, RHS.get(), G, *LHS, Y);
     X = *X1;

   }

   return status;
}

// This function solves
//    | J   A ||X|   |0|
//    | B^T C ||Y| = |G|
// via:  Xt = J^-1*A, Y = (C-B^T*Xt)^-1*G, X = -Xt*Y, where special cases of 
// C,G=0 are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Bordering::solveFZero(
		       Teuchos::ParameterList& params,
		       const NOX::Abstract::MultiVector* AA,
		       const LOCA::MultiContinuation::ConstraintInterface* BB,
		       const NOX::Abstract::MultiVector::DenseMatrix* CC,
		       const NOX::Abstract::MultiVector::DenseMatrix* G,
		       NOX::Abstract::MultiVector& X,
		       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSolver::Bordering::solveFZero()";
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
LOCA::BorderedSolver::Bordering::solveContiguous(
		       Teuchos::ParameterList& params,
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
    "LOCA::BorderedSolver::Bordering::solveContiguous()";
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

// This function solves
//    | J^T B ||X|   |0|
//    | A^T C ||Y| = |G|
// via:  Xt = J^-T*B, Y = (C-A^T*Xt)^-1*G, X = -Xt*Y, where special cases of 
// C,G=0 are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Bordering::solveFZeroTrans(
		       Teuchos::ParameterList& params,
		       const NOX::Abstract::MultiVector* AA,
		       const NOX::Abstract::MultiVector* BB,
		       const NOX::Abstract::MultiVector::DenseMatrix* CC,
		       const NOX::Abstract::MultiVector::DenseMatrix* G,
		       NOX::Abstract::MultiVector& X,
		       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSolver::Bordering::solveFTransZero()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Set X and Y to zero if G is zero
  if (isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return finalStatus;
  }

  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> Xt = 
    BB->clone(NOX::ShapeCopy);

  // compute Xt = J^-T B
  status = ts_grp->applyJacobianTransposeInverseMultiVector(params, *BB, *Xt);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  
  // compute t2 = -A^T*Xt
  NOX::Abstract::MultiVector::DenseMatrix t(AA->numVectors(),
					    Xt->numVectors());
  Xt->multiply(-1.0, *AA, t);
    
  // compute C^T - A^T*Xt
  if (!isZeroC)
    for (int i=0; i<t.numRows(); i++)
      for (int j=0; j<t.numCols(); j++)
	t(i,j) += (*CC)(j,i);

  // compute Y = (C^T - A^T*Xt)^-1 * G
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
// It also assumes F and B are in a contiguous multivec, stored in F
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::Bordering::solveContiguousTrans(
		       Teuchos::ParameterList& params,
		       const NOX::Abstract::MultiVector* AA,
		       const NOX::Abstract::MultiVector* BB,
		       const NOX::Abstract::MultiVector::DenseMatrix* CC,
		       vector<int>& indexF,
		       vector<int>& indexB,
		       const NOX::Abstract::MultiVector* F,
		       const NOX::Abstract::MultiVector::DenseMatrix* G,
		       NOX::Abstract::MultiVector& X,
		       NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSolver::Bordering::solveContiguousTrans()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // compute [X1 X2] = J^-T [F B]
  status = ts_grp->applyJacobianTransposeInverseMultiVector(params, *F, X);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X1 = X.subView(indexF);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> X2 = X.subView(indexB);
  
  // compute t1 = -A^T*X1, for efficiency t1 is stored in Y
  X1->multiply(-1.0, *AA, Y);
      
  // compute t2 = -A^T*X2
  NOX::Abstract::MultiVector::DenseMatrix t2(AA->numVectors(),
					     X2->numVectors());
  X2->multiply(-1.0, *AA, t2);

  // compute G - A^T*X1
  if (!isZeroG)
    Y += *G;
    
  // compute C^ - A^T*X2
  if (!isZeroC)
    for (int i=0; i<t2.numRows(); i++)
      for (int j=0; j<t2.numCols(); j++)
	t2(i,j) += (*CC)(j,i);

  // compute Y = (C^T - A^T*X2)^-1 * (G - A^T*X1)
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
