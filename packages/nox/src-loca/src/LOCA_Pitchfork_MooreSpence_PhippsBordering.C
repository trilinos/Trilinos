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

#include "LOCA_Pitchfork_MooreSpence_PhippsBordering.H"
#include "LOCA_Pitchfork_MooreSpence_ExtendedGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_LAPACK.hpp"  // for 4x4 matrix solve
#include "LOCA_BorderedSolver_JacobianOperator.H"

LOCA::Pitchfork::MooreSpence::PhippsBordering::PhippsBordering(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& slvrParams) : 
  globalData(global_data),
  solverParams(slvrParams),
  group(),
  pfGroup(),
  asymMultiVector(),
  asymVector(),
  nullVector(),
  JnVector(),
  dfdp(),
  dJndp(),
  borderedSolver(),
  nullMultiVector(),
  JnMultiVector(),
  sigma(0.0)
{
  borderedSolver = 
    globalData->locaFactory->createBorderedSolverStrategy(topParams,
							  solverParams);
}

LOCA::Pitchfork::MooreSpence::PhippsBordering::~PhippsBordering()
{
}

void
LOCA::Pitchfork::MooreSpence::PhippsBordering::setBlocks(
	 const Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::AbstractGroup>& group_,
	 const Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::ExtendedGroup>& pfGroup_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& asymMultiVector_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& nullVector_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& JnVector_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dfdp_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dJndp_)
{
  group = group_;
  pfGroup = pfGroup_;
  asymMultiVector = asymMultiVector_;
  asymVector = Teuchos::rcp(&(*asymMultiVector)[0], false);
  nullVector = nullVector_;
  JnVector = JnVector_;
  dfdp = dfdp_;
  dJndp = dJndp_;

  // Create multivectors for bordered solver
  nullMultiVector = nullVector->createMultiVector(1, NOX::DeepCopy);
  JnMultiVector = JnVector->createMultiVector(1, NOX::DeepCopy);
  sigma = JnVector->norm(NOX::Abstract::Vector::TwoNorm);
  JnMultiVector->scale(1.0/sigma);

  // Set blocks in bordered solver
  Teuchos::RefCountPtr<const LOCA::BorderedSolver::JacobianOperator> op =
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
LOCA::Pitchfork::MooreSpence::PhippsBordering::solve(
	   Teuchos::ParameterList& params,
	   const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& input,
           LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& result) const
{
  string callingFunction = 
    "LOCA::Pitchfork::MooreSpence::PhippsBordering::solve()";
  NOX::Abstract::Group::ReturnType status;
  
  // Get components of input
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> input_x = 
    input.getXMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> input_null = 
    input.getNullMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> input_slack = input.getSlacks();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> input_param = input.getBifParams();

  // Get components of result
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_x = 
    result.getXMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_null = 
    result.getNullMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> result_slack = 
    result.getSlacks();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    result.getBifParams();

  int m = input.numVectors();
  vector<int> index_input(m);
  for (int i=0; i<m; i++)
    index_input[i] = i;

  // Create new multivectors with m+3 columns
  // First m columns store input_x, input_null, result_x, result_null
  // respectively, next column stores dfdp, dJndp, J^-1 dfdp, J^-1 dJndp
  // respectively, next column stores psi, 0, J^-1 psi, etc...  
  // Last column is for solving (Jv)_x v
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+3);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_input_null = 
    input_null->clone(m+3);
  
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+3);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_null = 
    result_null->clone(m+3);

  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);
  
  // Set column m+1 to dfdp
  (*cont_input_x)[m] = *dfdp;
  
  // Set column m+2 to psi
  (*cont_input_x)[m+1] = *asymVector;
  
  // Initialize column m+3 to 0
  (*cont_input_x)[m+2].init(0.0);
  
  // Set first m columns to input_null
  cont_input_null->setBlock(*input_null, index_input);
  
  // Set column m+1 to dJndp
  (*cont_input_null)[m] = *dJndp;
  
  // Initialize column m+2 to 0
  (*cont_input_null)[m+1].init(0.0);
  
  // Initialize column m+3 to 0
  (*cont_input_null)[m+2].init(0.0);
  
  // Initialize result multivectors to 0
  cont_result_x->init(0.0);
  cont_result_null->init(0.0);
  
  // Solve
  status = solveContiguous(params, *cont_input_x, *cont_input_null, 
			   *input_slack, *input_param, 
			   *cont_result_x, *cont_result_null, 
			   *result_slack, *result_param);
  
  // Create views of first m columns for result_x, result_null
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_x_view = 
    cont_result_x->subView(index_input);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_null_view = 
    cont_result_null->subView(index_input);
  
  // Copy first m columns back into result_x, result_null
  *result_x = *cont_result_x_view;
  *result_null = *cont_result_null_view;

   return status;
}

// Solves pitchfork equations via Phipps modified bordering
// The first m columns of input_x and input_null store the RHS,
// column m+1 stores df/dp, d(Jn)/dp, column m+2 stores psi and 0,
// and the last column provides space for solving (Jv_x) v.  Note however
// input_param has only m columns.  result_x, result_null,
// are result_param have the same dimensions as their input counterparts
NOX::Abstract::Group::ReturnType 
LOCA::Pitchfork::MooreSpence::PhippsBordering::solveContiguous(
		  Teuchos::ParameterList& params,
		  const NOX::Abstract::MultiVector& input_x,
		  const NOX::Abstract::MultiVector& input_null,
		  const NOX::Abstract::MultiVector::DenseMatrix& input_slack,
	          const NOX::Abstract::MultiVector::DenseMatrix& input_param,
		  NOX::Abstract::MultiVector& result_x,
		  NOX::Abstract::MultiVector& result_null,
		  NOX::Abstract::MultiVector::DenseMatrix& result_slack,
	          NOX::Abstract::MultiVector::DenseMatrix& result_param) const
{
  string callingFunction = 
    "LOCA::Pitchfork::MooreSpence::PhippsBordering::solveContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  int m = input_x.numVectors()-3;
  vector<int> index_input(m);
  vector<int> index_input_dp(m+2);
  vector<int> index_null(1);
  vector<int> index_dp(1);
  vector<int> index_s(1);
  for (int i=0; i<m; i++) {
    index_input[i] = i;
    index_input_dp[i] = i;
  }
  index_input_dp[m] = m;
  index_input_dp[m+1] = m+1;
  index_dp[0] = m;
  index_s[0] = m+1;
  index_null[0] = m+2;

  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_1(1, m+2);
  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_2(1, m+3);

  // Create view of first m+2 columns of input_x, result_x
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> input_x_view = 
      input_x.subView(index_input_dp);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_x_view = 
      result_x.subView(index_input_dp);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  // Solve  |J   u||A B C| = |F df/dp psi|
  //        |v^T 0||a b c|   |0   0    0 |
  status = borderedSolver->applyInverse(params, input_x_view.get(), NULL, 
					*result_x_view, tmp_mat_1);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> A = 
    result_x.subView(index_input);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> B = 
    result_x.subView(index_dp);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> C = 
    result_x.subView(index_s);
  double b = tmp_mat_1(0,m);
  double c = tmp_mat_1(0,m+1);

  // compute (Jv)_x[A B C v]
  result_x[m+2] = *nullVector;
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> tmp = 
    result_x.clone(NOX::ShapeCopy);
  status = group->computeDJnDxaMulti(*nullVector, *JnVector, result_x,
				     *tmp);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // compute [G d(Jn)/dp 0 0] - (Jv)_x[A B C v]
  tmp->update(1.0, input_null, -1.0);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Solve  |J   u||D E K L| = |G-(Jv)_xA  d(Jv)/dp-(Jv)_xB  -(Jv)_xC -(Jv)_xv|
  //        |v^T 0||d e k l|   |    0             0              0        0   |
  status = borderedSolver->applyInverse(params, tmp.get(), NULL, result_null,
					tmp_mat_2);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> D = 
    result_null.subView(index_input);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> E = 
    result_null.subView(index_dp);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> K = 
    result_null.subView(index_s);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> L = 
    result_null.subView(index_null);
  double e = tmp_mat_2(0, m);
  double k = tmp_mat_2(0, m+1);
  double l = tmp_mat_2(0, m+2);

  double ltE = pfGroup->lTransNorm((*E)[0]);
  double ltK = pfGroup->lTransNorm((*K)[0]);
  double ltL = pfGroup->lTransNorm((*L)[0]);
  double ltv = pfGroup->lTransNorm(*nullVector);
  double ipv = group->innerProduct(*nullVector, *asymVector);
  double ipB = group->innerProduct((*B)[0], *asymVector);
  double ipC = group->innerProduct((*C)[0], *asymVector);

  // Fill coefficient arrays
  double M[16];
  M[0]  = sigma; M[1]  = -l;     M[2]  =  ipv; M[3]  =  ltL;
  M[4]  = 0.0;   M[5]  =  sigma; M[6]  =  0.0; M[7]  =  ltv;
  M[8]  = b;     M[9]  =  e;     M[10] = -ipB; M[11] = -ltE;
  M[12] = c;     M[13] =  k;     M[14] = -ipC; M[15] = -ltK;

  // compute s - <A,psi>
  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_3(1, m);
  group->innerProduct(*asymMultiVector, *A, tmp_mat_3);
  tmp_mat_3 -= input_slack;
  tmp_mat_3.scale(-1.0);

  // compute h - phi^T D
  NOX::Abstract::MultiVector::DenseMatrix tmp_mat_4(1, m);
  pfGroup->lTransNorm(*D, tmp_mat_4);
  tmp_mat_4 -= input_param;
  tmp_mat_4.scale(-1.0);

  double *R = new double[4*m];
  for (int i=0; i<m; i++) {
    R[4*i]   = tmp_mat_1(0,i);
    R[4*i+1] = tmp_mat_2(0,i);
    R[4*i+2] = tmp_mat_3(0,i);
    R[4*i+3] = tmp_mat_4(0,i);
  }

  // Solve M*P = R
  int piv[4];
  int info;
  Teuchos::LAPACK<int,double> dlapack;
  dlapack.GESV(4, m, M, 4, piv, R, 4, &info);
  if (info != 0) {
    globalData->locaErrorCheck->throwError(
				    callingFunction,
				    "Solve of 4x4 coefficient matrix failed!");
    return NOX::Abstract::Group::Failed;
  }

  NOX::Abstract::MultiVector::DenseMatrix alpha(1,m);
  NOX::Abstract::MultiVector::DenseMatrix beta(1,m);
  for (int i=0; i<m; i++) {
    alpha(0,i)        = R[4*i];
    beta(0,i)         = R[4*i+1];
    result_param(0,i) = R[4*i+2];
    result_slack(0,i) = R[4*i+3];
  }

  // compute A = A - B*z -C*w + v*alpha (remember A is a sub-view of result_x)
  A->update(Teuchos::NO_TRANS, -1.0, *B, result_param, 1.0);
  A->update(Teuchos::NO_TRANS, -1.0, *C, result_slack, 1.0);
  A->update(Teuchos::NO_TRANS, 1.0, *nullMultiVector, alpha, 1.0);

  // compute D = D - E*z - K*w + L*alpha + v*beta 
  // (remember D is a sub-view of result_null)
  D->update(Teuchos::NO_TRANS, -1.0, *E, result_param, 1.0);
  D->update(Teuchos::NO_TRANS, -1.0, *K, result_slack, 1.0);
  D->update(Teuchos::NO_TRANS, 1.0, *L, alpha, 1.0);
  D->update(Teuchos::NO_TRANS, 1.0, *nullMultiVector, beta, 1.0);

  delete [] R;

  return finalStatus;
}
