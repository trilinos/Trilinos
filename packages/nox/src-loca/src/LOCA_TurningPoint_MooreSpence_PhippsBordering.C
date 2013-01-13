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

#include "LOCA_TurningPoint_MooreSpence_PhippsBordering.H"
#include "LOCA_TurningPoint_MooreSpence_ExtendedGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_LAPACK.hpp"  // for 3x3 matrix solve
#include "LOCA_BorderedSolver_JacobianOperator.H"

LOCA::TurningPoint::MooreSpence::PhippsBordering::PhippsBordering(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& slvrParams) : 
  globalData(global_data),
  solverParams(slvrParams),
  group(),
  tpGroup(),
  nullVector(),
  JnVector(),
  dfdp(),
  dJndp(),
  borderedSolver(),
  nullMultiVector(),
  JnMultiVector(),
  s(0.0)
{
  borderedSolver = 
    globalData->locaFactory->createBorderedSolverStrategy(topParams,
							  solverParams);
}

LOCA::TurningPoint::MooreSpence::PhippsBordering::~PhippsBordering()
{
}

void
LOCA::TurningPoint::MooreSpence::PhippsBordering::setBlocks(
	 const Teuchos::RCP<LOCA::TurningPoint::MooreSpence::AbstractGroup>& group_,
	 const Teuchos::RCP<LOCA::TurningPoint::MooreSpence::ExtendedGroup>& tpGroup_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& nullVector_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& JnVector_,
	 const Teuchos::RCP<const NOX::Abstract::MultiVector>& dfdp_,
	 const Teuchos::RCP<const NOX::Abstract::MultiVector>& dJndp_)
{
  group = group_;
  tpGroup = tpGroup_;
  nullVector = nullVector_;
  JnVector = JnVector_;
  dfdp = dfdp_;
  dJndp = dJndp_;

  // Create multivectors for bordered solver
  nullMultiVector = nullVector->createMultiVector(1, NOX::DeepCopy);
  JnMultiVector = JnVector->createMultiVector(1, NOX::DeepCopy);
  s = JnVector->norm(NOX::Abstract::Vector::TwoNorm);
  JnMultiVector->scale(1.0/s);

  // Set blocks in bordered solver
  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> op =
    Teuchos::rcp(new  LOCA::BorderedSolver::JacobianOperator(group));
  borderedSolver->setMatrixBlocksMultiVecConstraint(op, 
						    JnMultiVector, 
						    nullMultiVector, 
						    Teuchos::null);
  NOX::Abstract::Group::ReturnType status = borderedSolver->initForSolve();
  globalData->locaErrorCheck->checkReturnType(status, 
		 "LOCA::Pitchfork::MooreSpence::PhippsBordering::setBlocks()");
}

NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::PhippsBordering::solve(
	   Teuchos::ParameterList& params,
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& input,
           LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::PhippsBordering::solve()";
  NOX::Abstract::Group::ReturnType status;
  
  // Get components of input
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x = 
    input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_null = 
    input.getNullMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = input.getScalars();

  // Get components of result
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x = 
    result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_null = 
    result.getNullMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    result.getScalars();

  int m = input.numVectors();
  std::vector<int> index_input(m);
  for (int i=0; i<m; i++)
    index_input[i] = i;
  
  // Create new multivectors with m+2 columns
  // First m columns store input_x, input_null, result_x, result_null
  // respectively, next column stores dfdp, dJndp, J^-1 dfdp, J^-1 dJndp
  // respectively.  Last column is for solving (Jv)_x v
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+2);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_null = 
    input_null->clone(m+2);
  
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+2);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_null = 
    result_null->clone(m+2);
  
  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);
  
  // Set column m+1 to dfdp
  (*cont_input_x)[m] = (*dfdp)[0];
  
  // Initialize column m+2 to 0
  (*cont_input_x)[m+1].init(0.0);
  
  // Set first m columns to input_null
  cont_input_null->setBlock(*input_null, index_input);
  
  // Set column m+1 to dJndp
  (*cont_input_null)[m] = (*dJndp)[0];
  
  // Initialize column m+2 to 0
  (*cont_input_null)[m+1].init(0.0);
  
  // Initialize result multivectors to 0
  cont_result_x->init(0.0);
  cont_result_null->init(0.0);
  
  // Solve
  status = solveContiguous(params, *cont_input_x, *cont_input_null, 
			   *input_param, *cont_result_x, *cont_result_null, 
			   *result_param);
  
  // Create views of first m columns for result_x, result_null
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x_view = 
    cont_result_x->subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_null_view = 
    cont_result_null->subView(index_input);
  
  // Copy first m columns back into result_x, result_null
  *result_x = *cont_result_x_view;
  *result_null = *cont_result_null_view;

  return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::PhippsBordering::solveTranspose(
	   Teuchos::ParameterList& params,
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& input,
           LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::PhippsBordering::solveTranspose()";
  NOX::Abstract::Group::ReturnType status;
  
  // Get components of input
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x = 
    input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_null = 
    input.getNullMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_param = input.getScalars();

  // Get components of result
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x = 
    result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_null = 
    result.getNullMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    result.getScalars();

  int m = input.numVectors();

  std::vector<int> index_input(m);
  for (int i=0; i<m; i++)
    index_input[i] = i;
  
  // Create new multivectors with m+2 columns
  // First m columns store input_x, input_null, result_x, result_null
  // respectively, next column stores 0, -phi, J^-T tmp , -J^-T phi
  // respectively, last column is for solving (Jv)_x^T u
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+2);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_null = 
    input_null->clone(m+2);
  
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+2);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_null = 
    result_null->clone(m+2);
  
  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);

  // Set last two columns to 0
  (*cont_input_x)[m].init(0);
  (*cont_input_x)[m+1].init(0);
  
    // Set first m columns to input_null
  cont_input_null->setBlock(*input_null, index_input);
  
  // Set next column to -phi
  Teuchos::RCP<NOX::Abstract::Vector> phi = tpGroup->getLengthVector();
  (*cont_input_null)[m].update(-1.0, *phi, 0.0);
  
  // Set last column to 0
  (*cont_input_null)[m].init(0);
  
  // Initialize result multivectors to 0
  cont_result_x->init(0.0);
  cont_result_null->init(0.0);
    
  // Solve
  status = solveTransposeContiguous(params, *cont_input_x, *cont_input_null, 
				    *input_param, *cont_result_x, 
				    *cont_result_null, 
				    *result_param);
  
  // Create views of first m columns for result_x, result_null
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x_view = 
    cont_result_x->subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_null_view = 
    cont_result_null->subView(index_input);
  
  // Copy first m columns back into result_x, result_null
  *result_x = *cont_result_x_view;
  *result_null = *cont_result_null_view;

   return status;
}

// Solves turning point equations via Phipps modified bordering
// The first m columns of input_x and input_null store the RHS while
// the last column stores df/dp, d(Jn)/dp respectively.  Note however
// input_param has only m columns (not m+1).  result_x, result_null,
// are result_param have the same dimensions as their input counterparts
NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::PhippsBordering::solveContiguous(
		  Teuchos::ParameterList& params,
		  const NOX::Abstract::MultiVector& input_x,
		  const NOX::Abstract::MultiVector& input_null,
	          const NOX::Abstract::MultiVector::DenseMatrix& input_param,
		  NOX::Abstract::MultiVector& result_x,
		  NOX::Abstract::MultiVector& result_null,
	          NOX::Abstract::MultiVector::DenseMatrix& result_param) const
{
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::PhippsBordering::solveContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  int m = input_x.numVectors()-2;
  std::vector<int> index_input(m);
  std::vector<int> index_input_dp(m+1);
  std::vector<int> index_null(1);
  std::vector<int> index_dp(1);
  for (int i=0; i<m; i++) {
    index_input[i] = i;
    index_input_dp[i] = i;
  }
  index_input_dp[m] = m;
  index_dp[0] = m;
  index_null[0] = m+1;

  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_1(1, m+1);
  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_2(1, m+2);

  // Create view of first m+1 columns of input_x, result_x
  Teuchos::RCP<NOX::Abstract::MultiVector> input_x_view = 
      input_x.subView(index_input_dp);
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x_view = 
      result_x.subView(index_input_dp);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  // Solve  |J   u||A B| = |F df/dp|
  //        |v^T 0||a b|   |0   0  |
  status = borderedSolver->applyInverse(params, input_x_view.get(), NULL, 
					*result_x_view, tmp_mat_1);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> A = 
    result_x.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> B = 
    result_x.subView(index_dp);
  double b = tmp_mat_1(0,m);

  // compute (Jv)_x[A B v]
  result_x[m+1] = *nullVector;
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp = 
    result_x.clone(NOX::ShapeCopy);
  status = group->computeDJnDxaMulti(*nullVector, *JnVector, result_x,
				     *tmp);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // compute (Jv)_x[A B v] - [G d(Jn)/dp 0]
  tmp->update(-1.0, input_null, 1.0);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Solve  |J   u||C D E| = |(Jv)_x A - G  (Jv)_x B - d(Jv)/dp  (Jv)_x v|
  //        |v^T 0||c d e|   |         0             0               0   |
  status = borderedSolver->applyInverse(params, tmp.get(), NULL, result_null,
					tmp_mat_2);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> C = 
    result_null.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> D = 
    result_null.subView(index_dp);
  Teuchos::RCP<NOX::Abstract::MultiVector> E = 
    result_null.subView(index_null);
  double d = tmp_mat_2(0, m);
  double e = tmp_mat_2(0, m+1);

  // Fill coefficient arrays
  double M[9];
  M[0] = s;   M[1] =  e;  M[2] = -tpGroup->lTransNorm((*E)[0]);
  M[3] = 0.0; M[4] =  s;  M[5] =  tpGroup->lTransNorm(*nullVector);
  M[6] = b;   M[7] = -d;  M[8] =  tpGroup->lTransNorm((*D)[0]);

  // compute h + phi^T C
  tpGroup->lTransNorm(*C, result_param);
  result_param += input_param;

  double *R = new double[3*m];
  for (int i=0; i<m; i++) {
    R[3*i]   =  tmp_mat_1(0,i);
    R[3*i+1] = -tmp_mat_2(0,i);
    R[3*i+2] =  result_param(0,i);
  }

  // Solve M*P = R
  int three = 3;
  int piv[3];
  int info;
  Teuchos::LAPACK<int,double> L;
  L.GESV(three, m, M, three, piv, R, three, &info);
  if (info != 0) {
    globalData->locaErrorCheck->throwError(
				    callingFunction,
				    "Solve of 3x3 coefficient matrix failed!");
    return NOX::Abstract::Group::Failed;
  }

  NOX::Abstract::MultiVector::DenseMatrix alpha(1,m);
  NOX::Abstract::MultiVector::DenseMatrix beta(1,m);
  for (int i=0; i<m; i++) {
    alpha(0,i)        = R[3*i];
    beta(0,i)         = R[3*i+1];
    result_param(0,i) = R[3*i+2];
  }

  // compute A = A - B*z + v*alpha (remember A is a sub-view of result_x)
  A->update(Teuchos::NO_TRANS, -1.0, *B, result_param, 1.0);
  A->update(Teuchos::NO_TRANS, 1.0, *nullMultiVector, alpha, 1.0);

  // compute C = -C + d*z - E*alpha + v*beta 
  // (remember C is a sub-view of result_null)
  C->update(Teuchos::NO_TRANS, 1.0, *D, result_param, -1.0);
  C->update(Teuchos::NO_TRANS, -1.0, *E, alpha, 1.0);
  C->update(Teuchos::NO_TRANS, 1.0, *nullMultiVector, beta, 1.0);

  delete [] R;

  return finalStatus;
}

// Solves turning point equations via classic Salinger bordering
// The first m columns of input_x and input_null store the RHS while
// the last column stores df/dp, d(Jn)/dp respectively.  Note however
// input_param has only m columns (not m+1).  result_x, result_null,
// are result_param have the same dimensions as their input counterparts
NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::PhippsBordering::solveTransposeContiguous(
		  Teuchos::ParameterList& params,
		  const NOX::Abstract::MultiVector& input_x,
		  const NOX::Abstract::MultiVector& input_null,
	          const NOX::Abstract::MultiVector::DenseMatrix& input_param,
		  NOX::Abstract::MultiVector& result_x,
		  NOX::Abstract::MultiVector& result_null,
	          NOX::Abstract::MultiVector::DenseMatrix& result_param) const
{
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::PhippsBordering::solveTransposeContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  int m = input_x.numVectors()-2;
  std::vector<int> index_input(m);
  std::vector<int> index_input_dp(m+1);
  std::vector<int> index_null(1);
  std::vector<int> index_dp(1);
  for (int i=0; i<m; i++) {
    index_input[i] = i;
    index_input_dp[i] = i;
  }
  index_input_dp[m] = m;
  index_dp[0] = m;
  index_null[0] = m+1;

  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_1(1, m+1);
  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_2(1, m+2);

  // Create view of first m+1 columns of input_null, result_null
  Teuchos::RCP<NOX::Abstract::MultiVector> input_null_view = 
      input_null.subView(index_input_dp);
  Teuchos::RCP<NOX::Abstract::MultiVector> result_null_view = 
      result_null.subView(index_input_dp);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Solve  |J^T v||A B| = |G -phi|
  //        |u^T 0||a b|   |0   0 |
  status =
    transposeBorderedSolver->applyInverseTranspose(params, 
						   input_null_view.get(), 
						   NULL, 
						   *result_null_view, 
						   tmp_mat_1);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> A = 
    result_null.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> B = 
    result_null.subView(index_dp);
  double b = tmp_mat_1(0,m);

  // compute (Jv)_x^T[A B u]
  result_null[m+1] = *uVector;
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp = 
    result_null.clone(NOX::ShapeCopy);
  status = group->computeDwtJnDxMulti(result_null, *nullVector, *tmp);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // compute [F 0 0] - (Jv)_x^T[A B u]
  tmp->update(1.0, input_x, -1.0);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Solve  |J^T v||C D E| = |F - (Jv)_x^T A  -(Jv)_x^T B  -(Jv)_x^T u|
  //        |u^T 0||c d e|   |         0             0            0   |
  status = 
    transposeBorderedSolver->applyInverseTranspose(params, 
						   tmp.get(), 
						   NULL, 
						   result_x,
						   tmp_mat_2);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> C = 
    result_x.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> D = 
    result_x.subView(index_dp);
  Teuchos::RCP<NOX::Abstract::MultiVector> E = 
    result_x.subView(index_null);
  double d = tmp_mat_2(0, m);
  double e = tmp_mat_2(0, m+1);

  // compute (Jv)_p^T*[A B u]
  NOX::Abstract::MultiVector::DenseMatrix t1(1,m+2);
  result_null.multiply(1.0, *dJndp, t1);

  // compute f_p^T*[C D E]
  NOX::Abstract::MultiVector::DenseMatrix t2(1,m+2);
  result_x.multiply(1.0, *dfdp, t2);

  // compute f_p^T*u
  double fptu = uVector->innerProduct((*dfdp)[0]);

  // Fill coefficient arrays
  double M[9];
  M[0] = st;   M[1] =  -e;   M[2] = t1(0,m+1) + t2(0,m+1);
  M[3] = 0.0;  M[4] =   st;  M[5] = fptu;
  M[6] = -b;   M[7] =  -d;   M[8] = t1(0,m) + t2(0,m);

  // Compute RHS
  double *R = new double[3*m];
  for (int i=0; i<m; i++) {
    R[3*i]   = tmp_mat_1(0,i);
    R[3*i+1] = tmp_mat_2(0,i);
    R[3*i+2] = result_param(0,i) - t1(0,i) - t2(0,i);
  }

  // Solve M*P = R
  int three = 3;
  int piv[3];
  int info;
  Teuchos::LAPACK<int,double> L;
  L.GESV(three, m, M, three, piv, R, three, &info);
  if (info != 0) {
    globalData->locaErrorCheck->throwError(
				    callingFunction,
				    "Solve of 3x3 coefficient matrix failed!");
    return NOX::Abstract::Group::Failed;
  }

  NOX::Abstract::MultiVector::DenseMatrix alpha(1,m);
  NOX::Abstract::MultiVector::DenseMatrix beta(1,m);
  for (int i=0; i<m; i++) {
    alpha(0,i)        = R[3*i];
    beta(0,i)         = R[3*i+1];
    result_param(0,i) = R[3*i+2];
  }

  // compute A = A + B*z + alpha*u (remember A is a sub-view of result_null)
  A->update(Teuchos::NO_TRANS, 1.0, *B, result_param, 1.0);
  A->update(Teuchos::NO_TRANS, 1.0, *uMultiVector, alpha, 1.0);

  // compute C = C + D*z + alpha*E + beta*u 
  // (remember C is a sub-view of result_x)
  C->update(Teuchos::NO_TRANS, 1.0, *D, result_param, 1.0);
  C->update(Teuchos::NO_TRANS, 1.0, *E, alpha, 1.0);
  C->update(Teuchos::NO_TRANS, 1.0, *uMultiVector, beta, 1.0);

  delete [] R;

  return finalStatus;
}
