// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include "LOCA_BorderedSolver_JacobianOperator.H"
#include "LOCA_BorderedSolver_ComplexOperator.H"
#include "LOCA_Hopf_ComplexMultiVector.H"

#ifdef HAVE_NOX_EPETRAEXT
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockMultiVector.h"
#endif

LOCA::BorderedSolver::EpetraHouseholder::EpetraHouseholder(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& slvrParams): 
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  op(),
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
  Ablock(),
  Bblock(),
  Ascaled(),
  Bscaled(),
  Cscaled(),
  linSys(),
  epetraOp(),
  baseMap(),
  globalMap(),
  numConstraints(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isValidForSolve(false),
  isValidForTransposeSolve(false),
  dblas(),
  scale_rows(true),
  scale_vals(),
  precMethod(JACOBIAN),
  includeUV(false),
  use_P_For_Prec(false),
  isComplex(false),
  omega(0.0)
{
  scale_rows = solverParams->get("Scale Augmented Rows", true);
  std::string prec_method = 
    solverParams->get("Preconditioner Method", "Jacobian");
  if (prec_method == "Jacobian")
    precMethod = JACOBIAN;
  else if (prec_method == "SMW")
    precMethod = SMW;
  else
    globalData->locaErrorCheck->throwError(
	    "LOCA::BorderedSolver::EpetraHouseholder::EpetraHouseholder()",
	    "Unknown preconditioner method!  Choices are Jacobian, SMW");
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
         const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op_,
	 const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
	 const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	 const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  std::string callingFunction = 
    "LOCA::BorderedSolver::EpetraHouseholder::setMatrixBlocks";

  op = op_;
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

    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tmpC = 
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

  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> jacOp =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::JacobianOperator>(op);
  Teuchos::RCP<const LOCA::BorderedSolver::ComplexOperator> complexOp =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::ComplexOperator>(op);

  if (jacOp != Teuchos::null) {

    isComplex = false;

    // Group must be an Epetra group
    Teuchos::RCP<const LOCA::Epetra::Group> constGrp = 
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(jacOp->getGroup());
    if (constGrp.get() == NULL)
      globalData->locaErrorCheck->throwError(
				    callingFunction,
				    "Group object must be an Epetra group");

    // Get A block
    Ablock = A;

    // Get B block
    Bblock = B;

    // Cast to non-const group
    grp = Teuchos::rcp_const_cast<LOCA::Epetra::Group>(constGrp);

    // Get linear system, and jacobian
    linSys = grp->getLinearSystem();
    epetraOp = linSys->getJacobianOperator();
  }
  else if (complexOp != Teuchos::null) {

    isComplex = true;

    omega = complexOp->getFrequency();

    // Group must be an Epetra group
    Teuchos::RCP<const LOCA::Epetra::Group> constGrp = 
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(complexOp->getGroup());
    if (constGrp.get() == NULL)
      globalData->locaErrorCheck->throwError(
				    callingFunction,
				    "Group object must be an Epetra group");

    // Cast to non-const group
    grp = Teuchos::rcp_const_cast<LOCA::Epetra::Group>(constGrp);

    // Get linear system and complex
    linSys = grp->getComplexLinearSystem();
    epetraOp = linSys->getJacobianOperator();
    
    // Get maps for complex vectors
    grp->getComplexMaps(baseMap, globalMap);

    // Get A block
    if (!isZeroA)
      Ablock = createBlockMV(*A);

    // Get A block
    if (!isZeroB)
      Bblock = createBlockMV(*B);
  }
  else {
    globalData->locaErrorCheck->throwError(
		      callingFunction,
		      std::string("Op argument must be of type !\n") + 
	              std::string("LOCA::BorderedSolver::JacobianOperator or \n") +
		      std::string("LOCA::BorderedSolver::ComplexOperator."));
  }

  Ascaled = Teuchos::null;
  Bscaled = Teuchos::null;
  Cscaled = Teuchos::null;
  
  if (Ablock != Teuchos::null)
    Ascaled = Ablock->clone();
  if (Bblock != Teuchos::null)
    Bscaled = Bblock->clone();
  if (C != Teuchos::null) {
    Cscaled = 
	  Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
								C->numRows(),
								C->numCols()));
    Cscaled->assign(*C);
  }
  
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

    // Scale rows
    if (scale_rows) {
      scale_vals.resize(numConstraints);
      double sn = std::sqrt( static_cast<double>(Bblock->length()) );
      for (int i=0; i<numConstraints; i++) {
	scale_vals[i] = (*Bblock)[i].norm() / sn;
	double t = 0.0;
	for (int j=0; j<numConstraints; j++)
	  t += (*C)(i,j) * (*C)(i,j);
	scale_vals[i] += std::sqrt(t);
	scale_vals[i] = 1.0 / scale_vals[i];
	(*Bscaled)[i].scale(scale_vals[i]);
	for (int j=0; j<numConstraints; j++)
	  (*Cscaled)(i,j) *= scale_vals[i];
      }
    }

    // Allocate vectors and matrices for factorization
    if (house_x.get() == NULL || house_x->numVectors() != numConstraints) {
      house_x = Bblock->clone(NOX::ShapeCopy);
      U = house_x->clone(NOX::ShapeCopy);
      V = house_x->clone(NOX::ShapeCopy);
      house_p.reshape(numConstraints, numConstraints);
      T.reshape(numConstraints, numConstraints);
      R.reshape(numConstraints, numConstraints);
    }

    // Factor constraints
    qrFact.computeQR(*Cscaled, *Bscaled, true, house_p, *house_x, T, R);

    // Compute U & V in P operator
    res = computeUV(house_p, *house_x, T, *Ablock, *U, *V, false);
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

    // Scale rows
    if (scale_rows) {
      scale_vals.resize(numConstraints);
      double sn = std::sqrt( static_cast<double>(Ablock->length()) );
      for (int i=0; i<numConstraints; i++) {
	scale_vals[i] = (*Ablock)[i].norm() / sn;
	double t = 0.0;
	for (int j=0; j<numConstraints; j++)
	  t += (*C)(j,i) * (*C)(j,i);
	scale_vals[i] += std::sqrt(t);
	scale_vals[i] = 1.0 / scale_vals[i];
	(*Ascaled)[i].scale(scale_vals[i]);
	for (int j=0; j<numConstraints; j++)
	  (*Cscaled)(j,i) *= scale_vals[i];
      }
    }

    // Allocate vectors and matrices for factorization
    if (house_x_trans.get() == NULL || 
	house_x_trans->numVectors() != numConstraints) {
      house_x_trans = Ablock->clone(NOX::ShapeCopy);
      U_trans = house_x_trans->clone(NOX::ShapeCopy);
      V_trans = house_x_trans->clone(NOX::ShapeCopy);
      house_p_trans.reshape(numConstraints, numConstraints);
      T_trans.reshape(numConstraints, numConstraints);
      R_trans.reshape(numConstraints, numConstraints);
    }

    // Factor constraints for transposed system
    qrFact.computeQR(*Cscaled, *Ascaled, false, house_p_trans, *house_x_trans, 
		     T_trans, R_trans);

    // Compute U & V in transposed P operator
    res = computeUV(house_p_trans, *house_x_trans, T_trans, *Bblock, *U_trans, 
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
  NOX::Abstract::Group::ReturnType status = op->apply(X, U);

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
  NOX::Abstract::Group::ReturnType status = op->applyTranspose(X, U);

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
    return ltbe.solve(params, *op, *constraints, *C, F, G, X, Y);
  }

  if (isZeroB) {
    LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
    return utbe.solve(params, *op, A.get(), *C, F, G, X, Y);
  }

  // Scale G
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> Gscaled; 
  if (G != NULL) {
    Gscaled = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
							       G->numCols()));
    if (scale_rows)
      for (int i=0; i<G->numRows(); i++)
	for (int j=0; j<G->numCols(); j++)
	  (*Gscaled)(i,j) = scale_vals[i]*(*G)(i,j);
    else
      Gscaled->assign(*G);
  }

  NOX::Abstract::Group::ReturnType res;
  if (!isComplex)
    res = solve(params, F, Gscaled.get(), X, Y);
  else {
    Teuchos::RCP<NOX::Abstract::MultiVector> blockF;
    if (!isZeroF)
      blockF = createBlockMV(*F);
    Teuchos::RCP<NOX::Abstract::MultiVector> blockX = createBlockMV(X);
    res = solve(params, blockF.get(), Gscaled.get(), *blockX, Y);
    setBlockMV(*blockX, X);
  }

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::EpetraHouseholder::applyInverseTranspose(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // If F and G are zero, return zero
  if (isZeroF && isZeroG) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  // If A or B is zero, we use bordering, which requires a transpose solve
  Teuchos::RCP<const LOCA::Abstract::TransposeSolveGroup> ts_grp;
  if (isZeroA || isZeroB) {
    ts_grp = 
      Teuchos::rcp_dynamic_cast<const LOCA::Abstract::TransposeSolveGroup>(grp);
    if (ts_grp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
	"LOCA::BorderedSolver::EpetraHouseholder::applyInverseTranspose()",
	string("Group must implement the LOCA::Abstract::TransposeSolveGroup")
	+ std::string(" interface in order to solve the transpose of the bordered")
	+ std::string(" system via bordering."));
  }

  if (isZeroA) { 
    LOCA::BorderedSolver::UpperTriangularBlockElimination utbe(globalData);
    return utbe.solveTranspose(params, *op, B.get(), *C, F, G, X, Y);
  }

  else if (isZeroB) { 
    LOCA::BorderedSolver::LowerTriangularBlockElimination ltbe(globalData);
    return ltbe.solveTranspose(params, *op, *A, *C, F, G, X, Y);
  }

  // Scale G
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> Gscaled; 
  if (G != NULL) {
    Gscaled = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
							       G->numCols()));
    if (scale_rows)
      for (int i=0; i<G->numRows(); i++)
	for (int j=0; j<G->numCols(); j++)
	  (*Gscaled)(i,j) = scale_vals[i]*(*G)(i,j);
    else
      Gscaled->assign(*G);
  }

  NOX::Abstract::Group::ReturnType res;
  if (!isComplex)
    res = solveTranspose(params, F, Gscaled.get(), X, Y);
  else {
    Teuchos::RCP<NOX::Abstract::MultiVector> blockF;
    if (!isZeroF)
      blockF = createBlockMV(*F);
    Teuchos::RCP<NOX::Abstract::MultiVector> blockX = createBlockMV(X);
    res = solveTranspose(params, blockF.get(), Gscaled.get(), *blockX, Y);
    setBlockMV(*blockX, X);
  }

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::EpetraHouseholder::solve(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  std::string callingFunction = 
    "LOCA::BorderedSolver::EpetraHouseholder::applyInverse()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // Make sure constraint factorization has been computed
  if (!isValidForSolve)
    globalData->locaErrorCheck->throwError(
		   callingFunction,
		   std::string("applyInverse() called with invalid constraint") + 
		   std::string(" factorizations.  Call initForSolve() first."));

  Teuchos::RCP<const NOX::Abstract::MultiVector> cRHS;

  if (!isZeroG) {
     
    Teuchos::RCP<NOX::Epetra::MultiVector> tmp_x = 
      Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(X.clone(NOX::ShapeCopy));
    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tmp_y = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
							       G->numCols()));
    Teuchos::RCP<NOX::Epetra::MultiVector> RHS = 
      Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(X.clone(NOX::ShapeCopy));

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
    int res = epetraOp->Apply(tmp_x->getEpetraMultiVector(), 
			      RHS->getEpetraMultiVector());
    if (res == 0)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::Failed;
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    RHS->update(Teuchos::NO_TRANS, -1.0, *Ablock, *tmp_y, -1.0);
    
    // Compute F - [A J]*P*[Z_y 0]
    if (!isZeroF) 
      RHS->update(1.0, *F, 1.0);
    
    cRHS = RHS;
  }
  else
    cRHS = Teuchos::rcp(F, false);
  
  // Get solution vec
  const NOX::Epetra::Vector& solution_vec = 
    dynamic_cast<const NOX::Epetra::Vector&>(grp->getX());
  
  // Create operator for P = J + U*V^T
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_U = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(U);
  Teuchos::RCP<Epetra_MultiVector> epetra_U = 
    Teuchos::rcp(&(nox_epetra_U->getEpetraMultiVector()), false);
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_V = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(V);
  Teuchos::RCP<Epetra_MultiVector> epetra_V = 
    Teuchos::rcp(&(nox_epetra_V->getEpetraMultiVector()), false);

  // Create a row-matrix version of P if J is a row matrix
  Teuchos::RCP<Epetra_RowMatrix> jac_rowmatrix = 
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetraOp);
  Teuchos::RCP<Epetra_Operator> op;
  if (jac_rowmatrix != Teuchos::null) {
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateRowMatrix(globalData, 
							       jac_rowmatrix, 
							       epetra_U, 
							       epetra_V,
							       false,
							       includeUV));
  }
  else
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData, 
							epetraOp, 
							epetra_U, 
							epetra_V,
							false));
  
  // Overwrite J with J + U*V^T if it's a CRS matrix and we aren't
  // using P for the preconditioner
  Teuchos::RCP<Epetra_CrsMatrix> jac_crs;
  if (includeUV && !use_P_For_Prec) {
    jac_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(epetraOp);
    if (jac_crs != Teuchos::null) {
      updateJacobianForPreconditioner(*U, *V, *jac_crs);
    }
  }

  // Set operator in solver to compute preconditioner
  if (use_P_For_Prec)
    linSys->setJacobianOperatorForSolve(op);
  else
    linSys->setJacobianOperatorForSolve(epetraOp);
     
  // Now compute preconditioner
  linSys->destroyPreconditioner();
  linSys->createPreconditioner(solution_vec, params, false);
   
  // Now recompute J if we modified it
  if (includeUV && !use_P_For_Prec && jac_crs != Teuchos::null) {
    grp->setX(solution_vec);
    if (isComplex)
      grp->computeComplex(omega);
    else
      grp->computeJacobian();
  }
       
  // Set operator for P in solver
  linSys->setJacobianOperatorForSolve(op);

  // Set preconditioner
  Teuchos::RCP<Epetra_Operator> prec_op;
  Teuchos::RCP<Epetra_Operator> epetraPrecOp;
  if (precMethod == SMW) {
    epetraPrecOp = linSys->getGeneratedPrecOperator();
    prec_op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData, 
							     epetraPrecOp, 
							     epetra_U, 
							     epetra_V,
							     true));
    linSys->setPrecOperatorForSolve(prec_op);
  }

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
  linSys->setJacobianOperatorForSolve(epetraOp);
  if (precMethod == SMW)
    linSys->setPrecOperatorForSolve(epetraPrecOp);
  linSys->destroyPreconditioner();

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSolver::EpetraHouseholder::solveTranspose(
			      Teuchos::ParameterList& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  std::string callingFunction = 
    "LOCA::BorderedSolver::EpetraHouseholder::applyInverseTranspose()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  bool isZeroF = (F == NULL);
  bool isZeroG = (G == NULL);

  // Make sure constraint factorization has been computed
  if (!isValidForTransposeSolve)
    globalData->locaErrorCheck->throwError(
	    callingFunction,
	    std::string("applyInverseTranspose() called with invalid constraint") + 
	    std::string(" factorizations.  Call initForSolve() first."));

  Teuchos::RCP<const NOX::Abstract::MultiVector> cRHS;
  
  if (!isZeroG) {
    
    Teuchos::RCP<NOX::Epetra::MultiVector> tmp_x = 
      Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(X.clone(NOX::ShapeCopy));
    Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tmp_y = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(G->numRows(),
							       G->numCols()));
    Teuchos::RCP<NOX::Epetra::MultiVector> RHS = 
      Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(X.clone(NOX::ShapeCopy));

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
    epetraOp->SetUseTranspose(true);
    int res = epetraOp->Apply(tmp_x->getEpetraMultiVector(), 
			      RHS->getEpetraMultiVector());
    epetraOp->SetUseTranspose(false);
    if (res == 0)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::Failed;
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    RHS->update(Teuchos::NO_TRANS, -1.0, *Bblock, *tmp_y, -1.0);

    // Compute F - [B J^T]*P*[Z_y 0]
    if (!isZeroF) 
      RHS->update(1.0, *F, 1.0);
    
    cRHS = RHS;
  }
  else
    cRHS = Teuchos::rcp(F, false);
  
  // Get solution vector
  const NOX::Epetra::Vector& solution_vec = 
    dynamic_cast<const NOX::Epetra::Vector&>(grp->getX());

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
//   Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy = tls_factory.create(Teuchos::rcp(&params, false), linSys);
  Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy = tls_factory.create(solverParams, linSys);
     
  // Compute Jacobian transpose (J^T)
  tls_strategy->createJacobianTranspose();
  Teuchos::RCP<Epetra_Operator> jac_trans = 
    tls_strategy->getJacobianTransposeOperator();

  // Create operator for P = J^T + U*V^T
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_U = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(U_trans);
  Teuchos::RCP<Epetra_MultiVector> epetra_U = 
    Teuchos::rcp(&(nox_epetra_U->getEpetraMultiVector()), false);
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_V = 
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(V_trans);
  Teuchos::RCP<Epetra_MultiVector> epetra_V = 
    Teuchos::rcp(&(nox_epetra_V->getEpetraMultiVector()), false);

  // Create a row-matrix version of P if J^T is a row matrix
  Teuchos::RCP<Epetra_RowMatrix> jac_trans_rowmatrix = 
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(jac_trans);
  Teuchos::RCP<Epetra_Operator> op;
  if (jac_trans_rowmatrix != Teuchos::null) 
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateRowMatrix(
							  globalData, 
							  jac_trans_rowmatrix, 
							  epetra_U, 
							  epetra_V,
							  false,
							  includeUV));
  else
    op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData, 
							jac_trans, 
							epetra_U, 
							epetra_V,
							false));
  
  // Overwrite J^T with J^T + U*V^T if it's a CRS matrix and we aren't
  // using P for the preconditioner
  Teuchos::RCP<Epetra_CrsMatrix> jac_trans_crs;
  if (includeUV && !use_P_For_Prec) {
    jac_trans_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(jac_trans);
    if (jac_trans_crs != Teuchos::null) {
      updateJacobianForPreconditioner(*U_trans, *V_trans, *jac_trans_crs);
    }
  }

   // Set operator in solver to compute preconditioner
  if (use_P_For_Prec)
    tls_strategy->setJacobianTransposeOperator(op);
  else
    tls_strategy->setJacobianTransposeOperator(jac_trans);
       
  // Now compute preconditioner
  tls_strategy->createTransposePreconditioner(solution_vec, params);
  
  // Now recompute J^T if we modified it
  if (includeUV && !use_P_For_Prec && jac_trans_crs != Teuchos::null) {
    grp->setX(solution_vec);
    if (isComplex)
      grp->computeComplex(omega);
    else
      grp->computeJacobian();
    tls_strategy->createJacobianTranspose();
  }
       
  // Set operator for P in transpose solver
  tls_strategy->setJacobianTransposeOperator(op);

  // Set preconditioner
  Teuchos::RCP<Epetra_Operator> prec_op;
  Teuchos::RCP<Epetra_Operator> epetraPrecOp;
  if (precMethod == SMW) {
    epetraPrecOp = tls_strategy->getTransposePreconditioner();
    prec_op = Teuchos::rcp(new LOCA::Epetra::LowRankUpdateOp(globalData, 
							     epetraPrecOp, 
							     epetra_U, 
							     epetra_V,
							     true));
    tls_strategy->setTransposePreconditioner(prec_op);
  }
  
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
       
  epetraOp->SetUseTranspose(false);

  // Set original operators in linear system
  linSys->setJacobianOperatorForSolve(epetraOp);
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
  const NOX::Epetra::MultiVector& Y2_epetra = 
    dynamic_cast<const NOX::Epetra::MultiVector&>(Y2);
  NOX::Epetra::MultiVector& UU_epetra = 
    dynamic_cast<NOX::Epetra::MultiVector&>(UU);
  bool use_trans = epetraOp->UseTranspose();
  epetraOp->SetUseTranspose(use_jac_transpose);
  int status = epetraOp->Apply(Y2_epetra.getEpetraMultiVector(), 
			       UU_epetra.getEpetraMultiVector());
  epetraOp->SetUseTranspose(use_trans);
  
  UU.update(Teuchos::NO_TRANS, 1.0, A, Y1, 1.0);
  VV.update(Teuchos::TRANS, 1.0, Y2, T, 0.0);

  if (status == 0)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::Failed;
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

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::BorderedSolver::EpetraHouseholder::createBlockMV(
				    const NOX::Abstract::MultiVector& v) const
{
#ifdef HAVE_NOX_EPETRAEXT
  const LOCA::Hopf::ComplexMultiVector& cv =
    dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(v);
  Teuchos::RCP<const NOX::Abstract::MultiVector> v_real =
    cv.getRealMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> v_imag =
    cv.getImagMultiVec();
  Teuchos::RCP<const NOX::Epetra::MultiVector> nox_epetra_v_real =
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(v_real);
  Teuchos::RCP<const NOX::Epetra::MultiVector> nox_epetra_v_imag =
    Teuchos::rcp_dynamic_cast<const NOX::Epetra::MultiVector>(v_imag);

  Teuchos::RCP<EpetraExt::BlockMultiVector> epetra_v = 
    Teuchos::rcp(new EpetraExt::BlockMultiVector(*baseMap, *globalMap,
						 v.numVectors()));
  epetra_v->LoadBlockValues(nox_epetra_v_real->getEpetraMultiVector(), 0);
  epetra_v->LoadBlockValues(nox_epetra_v_imag->getEpetraMultiVector(), 1);

  return Teuchos::rcp(new NOX::Epetra::MultiVector(
				      epetra_v,
				      NOX::DeepCopy,
				      NOX::Epetra::MultiVector::CreateView));
#else
  globalData->locaErrorCheck->throwError("LOCA::BorderedSolver::EpetraHouseholder::createBlockMV()", 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return Teuchos::null;
#endif
}

void
LOCA::BorderedSolver::EpetraHouseholder::setBlockMV(
				       const NOX::Abstract::MultiVector& bv,
				       NOX::Abstract::MultiVector& v) const
{
#ifdef HAVE_NOX_EPETRAEXT
  LOCA::Hopf::ComplexMultiVector& cv =
    dynamic_cast<LOCA::Hopf::ComplexMultiVector&>(v);
  Teuchos::RCP<NOX::Abstract::MultiVector> v_real =
    cv.getRealMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> v_imag =
    cv.getImagMultiVec();
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_v_real =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(v_real);
  Teuchos::RCP<NOX::Epetra::MultiVector> nox_epetra_v_imag =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(v_imag);

  const NOX::Epetra::MultiVector& nox_epetra_bv =
    dynamic_cast<const NOX::Epetra::MultiVector&>(bv);
  const Epetra_MultiVector& epetra_bv = nox_epetra_bv.getEpetraMultiVector();
  const EpetraExt::BlockMultiVector& block_bv =
    dynamic_cast<const EpetraExt::BlockMultiVector&>(epetra_bv);

  block_bv.ExtractBlockValues(nox_epetra_v_real->getEpetraMultiVector(), 0);
  block_bv.ExtractBlockValues(nox_epetra_v_imag->getEpetraMultiVector(), 1);
#else
  globalData->locaErrorCheck->throwError("LOCA::BorderedSolver::EpetraHouseholder::setBlockMV()", 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
#endif
}
