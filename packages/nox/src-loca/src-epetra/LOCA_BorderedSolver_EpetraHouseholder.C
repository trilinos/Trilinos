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

#include "LOCA_BorderedSolver_EpetraHouseholder.H"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "NOX_Epetra_MultiVector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_Epetra_Group.H"
#include "LOCA_Epetra_CompactWYOp.H"
#include "LOCA_Epetra_LowRankUpdateOp.H"
#include "LOCA_Epetra_LowRankUpdateRowMatrix.H"
#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H"
#include "LOCA_Epetra_TransposeLinearSystem_Factory.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_BorderedSolver_LowerTriangularBlockElimination.H"
#include "LOCA_BorderedSolver_UpperTriangularBlockElimination.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"

LOCA::BorderedSolver::EpetraHouseholder::EpetraHouseholder(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  A(),
  B(),
  C(),
  constraints(),
  qrFact(globalData),
  house_x(),
  house_p(),
  T(),
  R(),
  U(),
  V(),
  house_x_trans(),
  house_p_trans(),
  T_trans(),
  R_trans(),
  U_trans(),
  V_trans(),
  numConstraints(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isValidForSolve(false),
  isValidForTransposeSolve(false),
  dblas(),
  includeUV(false),
  use_P_For_Prec(false)
{
  includeUV = 
    solverParams->get("Include UV In Preconditioner", false);
  use_P_For_Prec = 
    solverParams->get("Use P For Preconditioner", false);
}

LOCA::BorderedSolver::EpetraHouseholder::~EpetraHouseholder()
{
}

void
LOCA::BorderedSolver::EpetraHouseholder::setMatrixBlocks(
         const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  string callingFunction = 
    "LOCA::BorderedSolver::EpetraHouseholder::setMatrixBlocks";

  // Group must be an Epetra group
  grp = Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(group);
  if (grp.get() == NULL)
    globalData->locaErrorCheck->throwError(
				    callingFunction,
				    "Group object must be an Epetra group");

  A = blockA;

  // Cast constraints to a ConstraintInterfaceMVDX
  constraints = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (constraints.get() == NULL)
    globalData->locaErrorCheck->throwError(
		 callingFunction,
		 "Constraints object must be of type ConstraintInterfaceMVDX");
  C = blockC;

  // Determine which blocks are zero
  isZeroA = (A.get() == NULL);
  isZeroB = constraints->isDXZero();
  isZeroC = (C.get() == NULL);

  // Get multivec constraint if it is nonzero
  if (isZeroB)
    B = Teuchos::null;
  else
    B = Teuchos::rcp(constraints->getDX(), false);

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

  // Set flags indicating we have to compute constraint factorizations
  isValidForSolve = false;
  isValidForTransposeSolve = false;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::EpetraHouseholder::initForSolve()
{
  NOX::Abstract::Group::ReturnType res = NOX::Abstract::Group::Ok;

  // Check if we need to compute anything
  if (isValidForSolve)
    return res;

  // Only compute QR factorization if A and B are nonzero
  if (!isZeroA && !isZeroB) {

    // Allocate vectors and matrices for factorization
    if (house_x.get() == NULL || house_x->numVectors() != numConstraints) {
      house_x = B->clone(NOX::ShapeCopy);
      U = house_x->clone(NOX::ShapeCopy);
      V = house_x->clone(NOX::ShapeCopy);
      house_p.reshape(numConstraints, numConstraints);
      T.reshape(numConstraints, numConstraints);
      R.reshape(numConstraints, numConstraints);
    }

    // Factor constraints
    qrFact.computeQR(*C, *B, true, house_p, *house_x, T, R);

    // Compute U & V in P operator
    res = computeUV(house_p, *house_x, T, *A, *U, *V, false);
    globalData->locaErrorCheck->checkReturnType(res, 
	       "LOCA::BorderedSolver::Epetra_Householder::initForSolve()");
  }

  isValidForSolve = true;

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::EpetraHouseholder::initForTransposeSolve()
{
  NOX::Abstract::Group::ReturnType res = NOX::Abstract::Group::Ok;

  // Check if we need to compute anything
  if (isValidForTransposeSolve)
    return res;

  // Only compute QR factorization if A and B are nonzero
  if (!isZeroA && !isZeroB) {

    // Allocate vectors and matrices for factorization
    if (house_x_trans.get() == NULL || 
	house_x_trans->numVectors() != numConstraints) {
      house_x_trans = A->clone(NOX::ShapeCopy);
      U_trans = house_x_trans->clone(NOX::ShapeCopy);
      V_trans = house_x_trans->clone(NOX::ShapeCopy);
      house_p_trans.reshape(numConstraints, numConstraints);
      T_trans.reshape(numConstraints, numConstraints);
      R_trans.reshape(numConstraints, numConstraints);
    }

    // Factor constraints for transposed system
    qrFact.computeQR(*C, *A, false, house_p_trans, *house_x_trans, T_trans, 
		     R_trans);

    // Compute U & V in transposed P operator
    res = computeUV(house_p_trans, *house_x_trans, T_trans, *B, *U_trans, 
		    *V_trans, true);
    globalData->locaErrorCheck->checkReturnType(res, 
	  "LOCA::BorderedSolver::Epetra_Householder::initForTransposeSolve()");
  }

  isValidForTransposeSolve = true;

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::EpetraHouseholder::apply(
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
LOCA::BorderedSolver::EpetraHouseholder::applyTranspose(
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
LOCA::BorderedSolver::EpetraHouseholder::applyInverse(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSolver::EpetraHouseholder::applyInverse()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F and G are zero, return zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  if (isZeroA) {
    LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
    return ltbe.solve(params, *grp, *constraints, *C, F, G, X, Y);
  }

  if (isZeroB) {
    LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
    return utbe.solve(params, *grp, A.get(), *C, F, G, X, Y);
  }

  // Make sure constraint factorization has been computed
  if (!isValidForSolve)
    globalData->locaErrorCheck->throwError(
		   callingFunction,
		   string("applyInverse() called with invalid constraint") + 
		   string(" factorizations.  Call initForSolve() first."));

  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> cRHS;

  if (!isZeroG) {
     
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> tmp_x = 
      X.clone(NOX::ShapeCopy);
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> tmp_y = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
							       G->numCols()));
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS = 
      X.clone(NOX::ShapeCopy);

    // Compute Z_y = R^-T * G
    Y.assign(*G);
    dblas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
	       Teuchos::NON_UNIT_DIAG,
	       G->numRows(), G->numCols(), 1.0, R.values(), 
	       R.numRows(), Y.values(), Y.numRows());
    
    // Compute P*[Z_y; 0]
    qrFact.applyCompactWY(house_p, *house_x, T, &Y, NULL, *tmp_y, *tmp_x, 
			  false);

    // Compute -[A J]*P*[Z_y 0]
    status = grp->applyJacobianMultiVector(*tmp_x, *RHS);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    RHS->update(Teuchos::NO_TRANS, -1.0, *A, *tmp_y, -1.0);
    
    // Compute F - [A J]*P*[Z_y 0]
    if (!isZeroF) 
      RHS->update(1.0, *F, 1.0);
    
    cRHS = RHS;
  }
  else
    cRHS = Teuchos::rcp(F, false);
   
  Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_a = 
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(A);
  Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_house_x = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(house_x);
  Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_a = 
    Teuchos::rcp(&(nox_epetra_a->getEpetraMultiVector()), false);
  Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_house_x = 
    Teuchos::rcp(&(nox_epetra_house_x->getEpetraMultiVector()), false);
  
  Teuchos::RefCountPtr<LOCA::Epetra::Group> nonconst_grp = 
    Teuchos::rcp_const_cast<LOCA::Epetra::Group>(grp);
     
  // Get linear system, jacobian, and solution vector
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> linSys = 
    nonconst_grp->getLinearSystem();
  Teuchos::RefCountPtr<Epetra_Operator> jac =
    linSys->getJacobianOperator();
  const NOX::Epetra::Vector& solution_vec = 
    dynamic_cast<const NOX::Epetra::Vector&>(grp->getX());
  
  // Create operator for P = J + U*V^T
  Teuchos::RefCountPtr<NOX::Epetra::MultiVector> nox_epetra_U = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(U);
  Teuchos::RefCountPtr<Epetra_MultiVector> epetra_U = 
    Teuchos::rcp(&(nox_epetra_U->getEpetraMultiVector()), false);
  Teuchos::RefCountPtr<NOX::Epetra::MultiVector> nox_epetra_V = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(V);
  Teuchos::RefCountPtr<Epetra_MultiVector> epetra_V = 
    Teuchos::rcp(&(nox_epetra_V->getEpetraMultiVector()), false);

  // Create a row-matrix version of P if J is a row matrix
  Teuchos::RefCountPtr<Epetra_RowMatrix> jac_rowmatrix = 
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(jac);
  Teuchos::RefCountPtr<Epetra_Operator> op;
  if (jac_rowmatrix != Teuchos::null) 
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateRowMatrix(globalData, 
							       jac_rowmatrix, 
							       epetra_U, 
							       epetra_V,
								includeUV));
  else
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData, 
							jac, 
							epetra_U, 
							epetra_V));
  
  // Overwrite J with J + U*V^T if it's a CRS matrix and we aren't
  // using P for the preconditioner
  Teuchos::RefCountPtr<Epetra_CrsMatrix> jac_crs;
  if (includeUV && !use_P_For_Prec) {
    jac_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(jac);
    if (jac_crs != Teuchos::null) {
      updateJacobianForPreconditioner(*U, *V, *jac_crs);
    }
  }

  // Set operator in solver to compute preconditioner
  if (use_P_For_Prec)
    linSys->setJacobianOperatorForSolve(op);
  else
    linSys->setJacobianOperatorForSolve(jac);
     
  // Now compute preconditioner
  linSys->destroyPreconditioner();
  linSys->createPreconditioner(solution_vec, params, false);
   
  // Now recompute J if we modified it
  if (includeUV && !use_P_For_Prec && jac_crs != Teuchos::null)
    linSys->computeJacobian(solution_vec);
       
  // Set operator for P in solver
  if (!use_P_For_Prec)
    linSys->setJacobianOperatorForSolve(op);
  
  // Solve for each RHS
  int m = X.numVectors();
  X.init(0.0);
  for (int i=0; i<m; i++) {
    bool stat = 
      linSys->applyJacobianInverse(
			 params, 
			 dynamic_cast<const NOX::Epetra::Vector&>((*cRHS)[i]),
			 dynamic_cast<NOX::Epetra::Vector&>(X[i]));
	 
    if (stat == true)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::NotConverged;
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  qrFact.applyCompactWY(house_p, *house_x, T, Y, X, isZeroG, false, false);

  // Set original operators in linear system
  linSys->setJacobianOperatorForSolve(jac);
  linSys->destroyPreconditioner();

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::EpetraHouseholder::applyInverseTranspose(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  string callingFunction = 
    "LOCA::BorderedSolver::EpetraHouseholder::applyInverseTranspose()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F and G are zero, return zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  // If A or B is zero, we use bordering, which requires a transpose solve
  Teuchos::RefCountPtr<const LOCA::Abstract::TransposeSolveGroup> ts_grp;
  if (isZeroA || isZeroB) {
    ts_grp = 
      Teuchos::rcp_dynamic_cast<const LOCA::Abstract::TransposeSolveGroup>(grp);
    if (ts_grp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
	 callingFunction,
	 string("Group must implement the LOCA::Abstract::TransposeSolveGroup")
	 + string(" interface in order to solve the transpose of the bordered")
	 + string(" system via bordering."));
  }

  if (isZeroA) { 
     LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
     return utbe.solveTranspose(params, *ts_grp, B.get(), *C, F, G, X, Y);
   }

   else if (isZeroB) { 
     LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
     return ltbe.solveTranspose(params, *ts_grp, *A, *C, F, G, X, Y);
   }

  // Make sure constraint factorization has been computed
  if (!isValidForTransposeSolve)
    globalData->locaErrorCheck->throwError(
	    callingFunction,
	    string("applyInverseTranspose() called with invalid constraint") + 
	    string(" factorizations.  Call initForSolve() first."));

  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> cRHS;
  
  if (!isZeroG) {
    
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> tmp_x = 
      X.clone(NOX::ShapeCopy);
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> tmp_y = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
							       G->numCols()));
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> RHS = 
      X.clone(NOX::ShapeCopy);

    // Compute Z_y = R^-T * G
    Y.assign(*G);
    dblas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
	       Teuchos::NON_UNIT_DIAG,
	       G->numRows(), G->numCols(), 1.0, R_trans.values(), 
	       R_trans.numRows(), Y.values(), Y.numRows());
    
    // Compute P*[Z_y; 0]
    qrFact.applyCompactWY(house_p_trans, *house_x_trans, T_trans, &Y, NULL, 
			  *tmp_y, *tmp_x, false);

    // Compute -[B J^T]*P*[Z_y 0]
    status = grp->applyJacobianTransposeMultiVector(*tmp_x, *RHS);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    RHS->update(Teuchos::NO_TRANS, -1.0, *B, *tmp_y, -1.0);

    // Compute F - [B J^T]*P*[Z_y 0]
    if (!isZeroF) 
      RHS->update(1.0, *F, 1.0);
    
    cRHS = RHS;
  }
  else
    cRHS = Teuchos::rcp(F, false);
   
  Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_b = 
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(B);
  Teuchos::RefCountPtr<const NOX::Epetra::MultiVector> nox_epetra_house_x = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(house_x_trans);
  Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_b = 
    Teuchos::rcp(&(nox_epetra_b->getEpetraMultiVector()), false);
  Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_house_x = 
    Teuchos::rcp(&(nox_epetra_house_x->getEpetraMultiVector()), false);
  
  Teuchos::RefCountPtr<LOCA::Epetra::Group> nonconst_grp = 
    Teuchos::rcp_const_cast<LOCA::Epetra::Group>(grp);
  
  // Get linear system, jacobian, and solution vector
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> linSys = 
    nonconst_grp->getLinearSystem();
  Teuchos::RefCountPtr<Epetra_Operator> jac =
    linSys->getJacobianOperator();
  const NOX::Epetra::Vector& solution_vec = 
    dynamic_cast<const NOX::Epetra::Vector&>(grp->getX());

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
//   Teuchos::RefCountPtr<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy = tls_factory.create(Teuchos::rcp(&params, false), linSys);
  Teuchos::RefCountPtr<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy = tls_factory.create(solverParams, linSys);
     
  // Compute Jacobian transpose (J^T)
  tls_strategy->computeJacobianTranspose(solution_vec);
  Teuchos::RefCountPtr<Epetra_Operator> jac_trans = 
    tls_strategy->getJacobianTransposeOperator();

  // Create operator for P = J^T + U*V^T
  Teuchos::RefCountPtr<NOX::Epetra::MultiVector> nox_epetra_U = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(U_trans);
  Teuchos::RefCountPtr<Epetra_MultiVector> epetra_U = 
    Teuchos::rcp(&(nox_epetra_U->getEpetraMultiVector()), false);
  Teuchos::RefCountPtr<NOX::Epetra::MultiVector> nox_epetra_V = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(V_trans);
  Teuchos::RefCountPtr<Epetra_MultiVector> epetra_V = 
    Teuchos::rcp(&(nox_epetra_V->getEpetraMultiVector()), false);

  // Create a row-matrix version of P if J^T is a row matrix
  Teuchos::RefCountPtr<Epetra_RowMatrix> jac_trans_rowmatrix = 
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(jac_trans);
  Teuchos::RefCountPtr<Epetra_Operator> op;
  if (jac_trans_rowmatrix != Teuchos::null) 
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateRowMatrix(
							  globalData, 
							  jac_trans_rowmatrix, 
							  epetra_U, 
							  epetra_V,
							  includeUV));
  else
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData, 
							jac_trans, 
							epetra_U, 
							epetra_V));
  
  // Overwrite J^T with J^T + U*V^T if it's a CRS matrix and we aren't
  // using P for the preconditioner
  Teuchos::RefCountPtr<Epetra_CrsMatrix> jac_trans_crs;
  if (includeUV && !use_P_For_Prec) {
    jac_trans_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(jac_trans);
    if (jac_trans_crs != Teuchos::null) {
      updateJacobianForPreconditioner(*U_trans, *V_trans, *jac_trans_crs);
    }
  }

   // Set operator in solver to compute preconditioner
  if (use_P_For_Prec)
    tls_strategy->setJacobianTransposeOperator(op);
       
  // Now compute preconditioner
  tls_strategy->createTransposePreconditioner(solution_vec, params);
  
  // Now recompute J^T if we modified it
  if (includeUV && !use_P_For_Prec && jac_trans_crs != Teuchos::null)
    tls_strategy->computeJacobianTranspose(solution_vec);
       
  // Set operator for P in transpose solver
  if (!use_P_For_Prec)
    tls_strategy->setJacobianTransposeOperator(op);
  
  // Solve for each RHS
  int m = X.numVectors();
  X.init(0.0);
  for (int i=0; i<m; i++) {
    bool stat = 
      tls_strategy->applyJacobianTransposeInverse(
			  params, 
			  dynamic_cast<const NOX::Epetra::Vector&>((*cRHS)[i]),
			  dynamic_cast<NOX::Epetra::Vector&>(X[i]));
    if (stat == true)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::NotConverged;
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

    
  }
	   
  qrFact.applyCompactWY(house_p_trans, *house_x_trans, T_trans, Y, X, isZeroG, 
			false, false);
       
  jac->SetUseTranspose(false);

  // Set original operators in linear system
  linSys->setJacobianOperatorForSolve(jac);
  linSys->destroyPreconditioner();

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraHouseholder::computeUV(
			    const NOX::Abstract::MultiVector::DenseMatrix& Y1,
			    const NOX::Abstract::MultiVector& Y2,
			    const NOX::Abstract::MultiVector::DenseMatrix& T,
			    const NOX::Abstract::MultiVector& A,
			    NOX::Abstract::MultiVector& UU,
			    NOX::Abstract::MultiVector& VV,
			    bool use_jac_transpose)
{
  // Now compute U and V in P = J + U*V^T, U = A*Y_1 + J*Y_2, V = Y_2*T^T
  NOX::Abstract::Group::ReturnType status;
  if (use_jac_transpose)
    status = grp->applyJacobianTransposeMultiVector(Y2, UU);
  else
    status = grp->applyJacobianMultiVector(Y2, UU);
  globalData->locaErrorCheck->checkReturnType(status, 
	       "LOCA::BorderedSolver::Epetra_Householder::computeUV()");
  
  UU.update(Teuchos::NO_TRANS, 1.0, A, Y1, 1.0);
  VV.update(Teuchos::TRANS, 1.0, Y2, T, 0.0);

  return status;
}

void
LOCA::BorderedSolver::EpetraHouseholder::updateJacobianForPreconditioner(
		      const NOX::Abstract::MultiVector& UU,
		      const NOX::Abstract::MultiVector& VV,
		      Epetra_CrsMatrix& jac) const
{
  // Get number of rows on this processor
  int numMyRows = jac.NumMyRows();

  int numEntries;
  int U_LDA;
  int V_LDA;
  int *row_indices;
  double *row_values;
  double *U_values;
  double *V_values;
  double val;

  // Get pointers to U & V values
  const NOX::Epetra::MultiVector nox_epetra_U = 
    dynamic_cast<const NOX::Epetra::MultiVector&>(UU);
  const Epetra_MultiVector& epetra_U = nox_epetra_U.getEpetraMultiVector();
  const NOX::Epetra::MultiVector& nox_epetra_V = 
    dynamic_cast<const NOX::Epetra::MultiVector&>(VV);
  const Epetra_MultiVector& epetra_V = nox_epetra_V.getEpetraMultiVector();
  epetra_U.ExtractView(&U_values, &U_LDA);
  epetra_V.ExtractView(&V_values, &V_LDA);

  const Epetra_BlockMap& v_map = epetra_V.Map();
  const Epetra_BlockMap& u_map = epetra_V.Map();
  const Epetra_BlockMap& row_map = jac.RowMap();
  const Epetra_BlockMap& col_map = jac.ColMap();

  for (int row=0; row<numMyRows; row++) {
    
    // extract view of row
    jac.ExtractMyRowView(row, numEntries, row_values, row_indices);

    int row_gid = row_map.GID(row);
    int u_row_lid = u_map.LID(row_gid);

    for (int col=0; col<numEntries; col++) {
      
      // Only included contributions from U*V^T on this proc
      int col_gid = col_map.GID(row_indices[col]);
      if (v_map.MyGID(col_gid)) {
	int v_row_lid = v_map.LID(col_gid);

	// val = sum_{k=0}^m U(i,k)*V(j,k)
	val = 0;
	for (int k=0; k<numConstraints; k++)
	  val += U_values[k*U_LDA+u_row_lid]*V_values[k*V_LDA+v_row_lid];

	// replace J(row,col) with J(row,col) + val
	row_values[col] += val;

      }
      
    }
  }
}
