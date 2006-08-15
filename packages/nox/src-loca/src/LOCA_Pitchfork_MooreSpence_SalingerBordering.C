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

#include "Teuchos_BLAS.hpp"
#include "LOCA_Pitchfork_MooreSpence_SalingerBordering.H"
#include "LOCA_Pitchfork_MooreSpence_ExtendedGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Pitchfork::MooreSpence::SalingerBordering::SalingerBordering(
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
  dJndp()
{
}

LOCA::Pitchfork::MooreSpence::SalingerBordering::~SalingerBordering()
{
}

void
LOCA::Pitchfork::MooreSpence::SalingerBordering::setBlocks(
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
}

NOX::Abstract::Group::ReturnType 
LOCA::Pitchfork::MooreSpence::SalingerBordering::solve(
	   Teuchos::ParameterList& params,
	   const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& input,
           LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& result) const
{
  string callingFunction = 
    "LOCA::Pitchfork::MooreSpence::SalingerBordering::solve()";
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

  // Create new multivectors with m+2 columns
  // First m columns store input_x, input_null, result_x, result_null
  // respectively, next column stores dfdp, dJndp, J^-1 dfdp, J^-1 dJndp
  // respectively, last column stores psi, 0, J^-1 psi, etc...
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+2);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_input_null = 
    input_null->clone(m+2);
  
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+2);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_null = 
    result_null->clone(m+2);
  
  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);
  
  // Set column m+1 to dfdp
  (*cont_input_x)[m] = *dfdp;
  
  // Set column m+2 to psi
  (*cont_input_x)[m+1] = *asymVector;
  
  // Set first m columns to input_null
  cont_input_null->setBlock(*input_null, index_input);
  
  // Set column m+1 to dJndp
  (*cont_input_null)[m] = *dJndp;
  
  // Set column m+2 to 0
  (*cont_input_null)[m+1].init(0.0);
  
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

// Solves pitchfork equations via classic Salinger bordering
// The first m columns of input_x and input_null store the RHS,
// column m+1 stores df/dp, d(Jn)/dp and column m+2 stores psi and 0 
// respectively.  Note however input_slack input_param have only m columns 
// (not m+2).  result_x, result_null, result_slack, and result_param have the 
// same dimensions as their input counterparts
NOX::Abstract::Group::ReturnType 
LOCA::Pitchfork::MooreSpence::SalingerBordering::solveContiguous(
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
    "LOCA::Pitchfork::MooreSpence::SalingerBordering::solveContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  int m = input_x.numVectors()-2;
  vector<int> index_input(m);
  vector<int> index_dp(1);
  vector<int> index_s(1);
  for (int i=0; i<m; i++)
    index_input[i] = i;
  index_dp[0] = m;
  index_s[0] = m+1;

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  // compute [A b c] = J^-1 [F df/dp psi]
  status = group->applyJacobianInverseMultiVector(params, input_x, result_x);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> A = 
    result_x.subView(index_input);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> b = 
    result_x.subView(index_dp);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> c = 
    result_x.subView(index_s);

  // compute (Jn)_x[A b c]
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> tmp = 
    result_x.clone(NOX::ShapeCopy);
  status = group->computeDJnDxaMulti(*nullVector, *JnVector, result_x,
				     *tmp);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // compute [G d(Jn)/dp 0] - (Jn)_x[A b c]
  tmp->update(1.0, input_null, -1.0);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute [D e g] = J^-1 [G d(Jn)/dp 0] - (Jn)_x[A b c]
  status = group->applyJacobianInverseMultiVector(params, *tmp, result_null);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> D = 
    result_null.subView(index_input);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> e = 
    result_null.subView(index_dp);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> g = 
    result_null.subView(index_s);

  // compute w = (phi^T e)(<A,psi> - s) - <b,psi>(phi^T D - h) / 
  //                ( (phi^T e) <c,psi> - (phi^T g) <b,psi> )
  Teuchos::BLAS<int,double> dblas;
  double lte = pfGroup->lTransNorm((*e)[0]);
  double ltf = pfGroup->lTransNorm((*g)[0]);
  double ipb = group->innerProduct((*b)[0], *asymVector);
  double ipc = group->innerProduct((*c)[0], *asymVector);
  double denom = lte*ipc - ltf*ipb;
  group->innerProduct(*asymMultiVector, *A, result_slack);
  pfGroup->lTransNorm(*D, result_param);
 
  for (int i=0; i<m; i++) {
    result_slack(0,i) = (lte*(result_slack(0,i) - input_slack(0,i)) - 
			 ipb*(result_param(0,i) - input_param(0,i))) / denom;
    result_param(0,i) = (result_param(0,i) - input_param(0,i) - 
			 ltf*result_slack(0,i)) / lte;
  }

  // compute A = A - b*z - c*w (remember A is a sub-view of result_x)
  A->update(Teuchos::NO_TRANS, -1.0, *b, result_param, 1.0);
  A->update(Teuchos::NO_TRANS, -1.0, *c, result_slack, 1.0);

  // compute D = D - e*z - g*w (remember D is a sub-view of result_null)
  D->update(Teuchos::NO_TRANS, -1.0, *e, result_param, 1.0);
  D->update(Teuchos::NO_TRANS, -1.0, *g, result_slack, 1.0);

  return finalStatus;
}
