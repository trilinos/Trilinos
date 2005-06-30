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

#include "LOCA_BorderedSystem_EpetraHouseholder.H"
#include "Epetra_MultiVector.h"
#include "NOX_Epetra_MultiVector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_EpetraNew_Group.H"
#include "LOCA_Epetra_CompactWYOp.H"

LOCA::BorderedSystem::EpetraHouseholder::EpetraHouseholder(
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
  house_x(),
  house_p(),
  T(),
  R(),
  v_x(),
  numConstraints(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isContiguous(false),
  dblas()
{
}

LOCA::BorderedSystem::EpetraHouseholder::~EpetraHouseholder()
{
}

void
LOCA::BorderedSystem::EpetraHouseholder::setIsContiguous(bool flag)
{
  isContiguous = flag;
}

void
LOCA::BorderedSystem::EpetraHouseholder::setMatrixBlocks(
	 const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  string callingFunction = 
    "LOCA::BorderedSystem::EpetraHouseholder::setMatrixBlocks";

  // Cast group to an EpetraNew group
  grp = Teuchos::rcp_dynamic_cast<const LOCA::EpetraNew::Group>(group);
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

  // We only use the Householder technique if A and B are nonzero 
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

  // factor constraints
  if (!isZeroA && !isZeroB)
    factorConstraints();

}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::EpetraHouseholder::apply(
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
LOCA::BorderedSystem::EpetraHouseholder::applyTranspose(
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
LOCA::BorderedSystem::EpetraHouseholder::applyInverse(
			      NOX::Parameter::List& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSystem::EpetraHouseholder::applyInverse()";
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

     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS;
     Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> cRHS;
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector> tmp_x;
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> Z_y;
     Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> tmp_y;

     if (!isZeroG) {

       tmp_x = Teuchos::rcp(x->clone(NOX::ShapeCopy));
       tmp_y = 
	 Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							       G->numRows(),
							       G->numCols()));
       RHS = Teuchos::rcp(x->clone(NOX::ShapeCopy));

       // Compute Z_y = R^-T * G
       Y.assign(*G);
       dblas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		  Teuchos::NON_UNIT_DIAG,
		  G->numRows(), G->numCols(), 1.0, R.values(), 
		  R.numRows(), Y.values(), Y.numRows());

       // Compute P*[Z_y; 0]
       applyCompactWY(NULL, &Y, *tmp_x, *tmp_y, false);

       // Compute -[A J]*P*[Z_y 0]
       status = grp->applyJacobianMultiVector(*tmp_x, *RHS);
       finalStatus = 
	 LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						      callingFunction);
       RHS->update(Teuchos::NO_TRANS, -1.0, *A, *tmp_y, -1.0);

       // Compute F - [A J]*P*[Z_y 0]
       if (!isZeroF) 
	 RHS->update(1.0, *cf, 1.0);

       cRHS = RHS;
     }
     else
       cRHS = cf;
   
     Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_a = 
       Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(ca);
     Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_x = 
       Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(x);
     Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_house_x = 
       Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(house_x);
     Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_a = 
       Teuchos::rcp(&(nox_epetra_a->getEpetraMultiVector()), false);
     Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_x = 
       Teuchos::rcp(&(nox_epetra_x->getEpetraMultiVector()), false);
      Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_house_x = 
       Teuchos::rcp(&(nox_epetra_house_x->getEpetraMultiVector()), false);
     
     const NOX::EpetraNew::LinearSystem& linSys = grp->getLinearSystem();
     Teuchos::RefCountPtr<const Epetra_Operator> jac =
       Teuchos::rcp(&linSys.getJacobianOperator(), false);

     LOCA::Epetra::CompactWYOp op(jac, epetra_a, epetra_house_x,
				  Teuchos::rcp(&house_p, false),
				  Teuchos::rcp(&T, false));

     op.init(*epetra_x);

     grp->setJacobianOperatorForSolve(op);

     status = grp->applyJacobianInverseMultiVector(params, *cRHS, *x);
     finalStatus = 
       LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						    callingFunction);

     applyCompactWY(*x, Y, false, isZeroG, false);

     op.finish();

     grp->setJacobianOperatorForSolve(*jac);
   }

   return finalStatus;

}

void
LOCA::BorderedSystem::EpetraHouseholder::factorConstraints()
{
  double beta;

  // Allocate house_x, house_p, beta if necesary, and copy dg/dx into house_x
  if (house_x.get() == NULL || house_x->numVectors() != numConstraints) {
    house_x = Teuchos::rcp(B->clone(NOX::DeepCopy));
    house_p.reshape(numConstraints, numConstraints);
    house_p.putScalar(0.0);
    T.reshape(numConstraints, numConstraints);
    T.putScalar(0.0);
    R.reshape(numConstraints, numConstraints);
    v_x = Teuchos::rcp(B->clone(1));
  }
  else 
    *house_x = *B;

  // Copy transpose of dg/dp into R
  for (int i=0; i<numConstraints; i++)
    for (int j=0; j<numConstraints; j++)
      R(i,j) = (*C)(j,i);

  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> v_p;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> h_x;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> h_p;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> y_x;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> y_p;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> z;
  vector<int> h_idx;
  vector<int> y_idx;
  y_idx.reserve(numConstraints);

  for (int i=0; i<numConstraints; i++) {

    // Create view of column i of house_p starting at row i
    v_p = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							 Teuchos::View, 
							 house_p, 
							 numConstraints-i, 
							 1, i, i));

    // Create view of columns i through numConstraints-1 of house_x
    h_idx.resize(numConstraints-i);
    for (unsigned int j=0; j<h_idx.size(); j++)
      h_idx[j] = i+j;
    h_x = Teuchos::rcp(house_x->subView(h_idx));

    // Create view of columns i thru numConstraints-1 of R, starting at row i
    h_p = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							 Teuchos::View, 
							 R,
							 numConstraints-i,
							 numConstraints-i,
							 i, i));

    if (i > 0) {

      // Create view of columns 0 through i-1 of house_x
      y_idx.push_back(i-1);
      y_x = Teuchos::rcp(house_x->subView(y_idx));
      
      // Create view of columns 0 through i-1 of house_p, starting at row i
      y_p = 
	Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							 Teuchos::View, 
							 house_p,
							 numConstraints-i,
							 i, i, 0));

      // Create view of column i, row 0 through i-1 of T
      z = 
	Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							 Teuchos::View, 
							 T, i, 1, 0, i));
    }

    // Compute Householder Vector
    computeHouseholderVector(i, *house_x, R, *v_x, *v_p, beta);

    // Apply Householder reflection
    applyHouseholderVector(*v_x, *v_p, beta, *h_x, *h_p);

    // Copy v_x into house_x
    (*house_x)[i] = (*v_x)[0];
    
    T(i,i) = -beta;

    if (i > 0) {

      // Compute z = y_x^T * v_x
      v_x->multiply(1.0, *y_x, *z);

      // Compute z = -beta * (z + y_p^T * v_p)
      z->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -beta, *y_p, *v_p, -beta);

      // Compute z = T * z
      dblas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
		 i, T.values(), numConstraints, z->values(), 1);

    }
  }

}

void
LOCA::BorderedSystem::EpetraHouseholder::computeHouseholderVector(
			  int col,
			  const NOX::Abstract::MultiVector& A_x,
			  const NOX::Abstract::MultiVector::DenseMatrix& A_p,
			  NOX::Abstract::MultiVector& V_x,
			  NOX::Abstract::MultiVector::DenseMatrix& V_p,
			  double& beta)
{
  double houseP = A_p(col,col);
  
  V_p(0,0) = 1.0;
  V_x[0] = A_x[col];

  double sigma = A_x[col].dot(A_x[col]);
  for (int i=col+1; i<A_p.numRows(); i++)    
    sigma += A_p(i,col)*A_p(i,col);

  if (sigma == 0.0)
    beta = 0.0;
  else {
    double mu = sqrt(houseP*houseP + sigma);
    if (houseP <= 0.0)
      houseP = houseP - mu;
    else
      houseP = -sigma / (houseP + mu);
    beta = 2.0*houseP*houseP/(sigma + houseP*houseP);
    
    V_x.scale(1.0/houseP);
    for (int i=1; i<V_p.numRows(); i++)
      V_p(i,0) = A_p(i+col,col) / houseP;
  }


  return;
}

void
LOCA::BorderedSystem::EpetraHouseholder::applyHouseholderVector(
			   const NOX::Abstract::MultiVector& V_x,
			   const NOX::Abstract::MultiVector::DenseMatrix& V_p,
			   double beta,
			   NOX::Abstract::MultiVector& A_x,
			   NOX::Abstract::MultiVector::DenseMatrix& A_p)
{
  int nColsA = A_x.numVectors();

  // Compute u = V_x^T * A_x
  NOX::Abstract::MultiVector::DenseMatrix u(1, nColsA);
  A_x.multiply(1.0, V_x, u);

  // Compute u = u + V_p^T * A_P
  u.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, V_p, A_p, 1.0);

  // Compute A_p = A_p - b*V_p*u
  A_p.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -beta, V_p, u, 1.0);

  // Compute A_x = A_x - b*V_x*u
  A_x.update(Teuchos::NO_TRANS, -beta, V_x, u, 1.0);
}

void
LOCA::BorderedSystem::EpetraHouseholder::applyCompactWY(
				   NOX::Abstract::MultiVector& x,
				   NOX::Abstract::MultiVector::DenseMatrix& y,
				   bool isZeroX, bool isZeroY,
				   bool useTranspose) const
{
  if (isZeroX && isZeroY) {
    x.init(0.0);
    y.putScalar(0.0);
    return;
  }

  Teuchos::ETransp T_flag;
  if (useTranspose)
    T_flag = Teuchos::TRANS;
  else
    T_flag = Teuchos::NO_TRANS;

  NOX::Abstract::MultiVector::DenseMatrix tmp(numConstraints, x.numVectors());

  // Compute Y_p^T*y + Y_x^T*x
  if (!isZeroX)
    x.multiply(1.0, *house_x, tmp);

  // Opportunity for optimization here since house_p is a lower-triangular
  // matrix with unit diagonal
  if (!isZeroX && !isZeroY)
    tmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, house_p, y, 1.0);
  else if (!isZeroY)
    tmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, house_p, y, 0.0);

  // Compute op(T)*(Y_p^T*y + Y_x^T*x)
  dblas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, T_flag, 
	     Teuchos::NON_UNIT_DIAG, tmp.numRows(), tmp.numCols(), 1.0, 
	     T.values(), T.numRows(), tmp.values(), tmp.numRows());

  // Compute y = y + Y_p*op(T)*(Y_p^T*y + Y_x^T*x)
  // Opportunity for optimization here since house_p is a lower-triangular
  // matrix with unit diagonal
  if (isZeroY)
    y.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, house_p, tmp, 0.0);
  else
    y.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, house_p, tmp, 1.0);

  // Compute x = x + Y_p*op(T)*(Y_p^T*y + Y_x^T*x)
  if (isZeroX)
    x.update(Teuchos::NO_TRANS, 1.0, *house_x, tmp, 0.0);
  else
    x.update(Teuchos::NO_TRANS, 1.0, *house_x, tmp, 1.0); 
}

void
LOCA::BorderedSystem::EpetraHouseholder::applyCompactWY(
		     const NOX::Abstract::MultiVector* input_x,
		     const NOX::Abstract::MultiVector::DenseMatrix* input_y,
		     NOX::Abstract::MultiVector& result_x,
		     NOX::Abstract::MultiVector::DenseMatrix& result_y,
		     bool useTranspose) const
{
  bool isZeroX = (input_x == NULL);
  bool isZeroY = (input_y == NULL);

  if (!isZeroX)
    result_x = *input_x;
  if (!isZeroY)
    result_y.assign(*input_y);

  applyCompactWY(result_x, result_y, isZeroX, isZeroY, useTranspose);
}

// This function solves
//    | J   0 ||X|   |F|
//    | B^T C ||Y| = |G|
// via:  X = J^-1 * F, Y = C^-1 * (G - B^T*X), where special cases of B,F,G=0
// are taken into account
NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::EpetraHouseholder::solveAZero(
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
LOCA::BorderedSystem::EpetraHouseholder::solveBZero(
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
