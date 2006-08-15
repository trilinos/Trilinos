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

#include "LOCA_TurningPoint_MooreSpence_SalingerBordering.H"
#include "LOCA_TurningPoint_MooreSpence_ExtendedGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::TurningPoint::MooreSpence::SalingerBordering::SalingerBordering(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& slvrParams) : 
  globalData(global_data),
  solverParams(slvrParams),
  group(),
  tpGroup(),
  nullVector(),
  JnVector(),
  dfdp(),
  dJndp()
{
}

LOCA::TurningPoint::MooreSpence::SalingerBordering::~SalingerBordering()
{
}

void
LOCA::TurningPoint::MooreSpence::SalingerBordering::setBlocks(
	 const Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::AbstractGroup>& group_,
	 const Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::ExtendedGroup>& tpGroup_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& nullVector_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& JnVector_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dfdp_,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dJndp_)
{
  group = group_;
  tpGroup = tpGroup_;
  nullVector = nullVector_;
  JnVector = JnVector_;
  dfdp = dfdp_;
  dJndp = dJndp_;
}

NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::SalingerBordering::solve(
	   Teuchos::ParameterList& params,
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& input,
           LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& result) const
{
  string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::SalingerBordering::solve()";
  NOX::Abstract::Group::ReturnType status;
  
  // Get components of input
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> input_x = 
    input.getXMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> input_null = 
    input.getNullMultiVec();
  Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> input_param = input.getScalars();

  // Get components of result
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_x = 
    result.getXMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> result_null = 
    result.getNullMultiVec();
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector::DenseMatrix> result_param = 
    result.getScalars();

  int m = input.numVectors();

  vector<int> index_input(m);
  for (int i=0; i<m; i++)
    index_input[i] = i;
  
  // Create new multivectors with m+1 columns
  // First m columns store input_x, input_null, result_x, result_null
  // respectively, last column stores dfdp, dJndp, J^-1 dfdp, J^-1 dJndp
  // respectively
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+1);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_input_null = 
    input_null->clone(m+1);
  
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+1);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> cont_result_null = 
    result_null->clone(m+1);
  
  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);

  // Set last column to dfdp
  (*cont_input_x)[m] = *dfdp;
  
    // Set first m columns to input_null
  cont_input_null->setBlock(*input_null, index_input);
  
  // Set last column to dJndp
  (*cont_input_null)[m] = *dJndp;
  
  // Initialize result multivectors to 0
  cont_result_x->init(0.0);
  cont_result_null->init(0.0);
    
  // Solve
  status = solveContiguous(params, *cont_input_x, *cont_input_null, 
			   *input_param, *cont_result_x, *cont_result_null, 
			   *result_param);
  
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

// Solves turning point equations via classic Salinger bordering
// The first m columns of input_x and input_null store the RHS while
// the last column stores df/dp, d(Jn)/dp respectively.  Note however
// input_param has only m columns (not m+1).  result_x, result_null,
// are result_param have the same dimensions as their input counterparts
NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::SalingerBordering::solveContiguous(
		  Teuchos::ParameterList& params,
		  const NOX::Abstract::MultiVector& input_x,
		  const NOX::Abstract::MultiVector& input_null,
	          const NOX::Abstract::MultiVector::DenseMatrix& input_param,
		  NOX::Abstract::MultiVector& result_x,
		  NOX::Abstract::MultiVector& result_null,
	          NOX::Abstract::MultiVector::DenseMatrix& result_param) const
{
  string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::SalingerBordering::solveContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  int m = input_x.numVectors()-1;
  vector<int> index_input(m);
  vector<int> index_dp(1);
  for (int i=0; i<m; i++)
    index_input[i] = i;
  index_dp[0] = m;

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  // compute [A b] = J^-1 [F df/dp]
  status = group->applyJacobianInverseMultiVector(params, input_x, result_x);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> A = 
    result_x.subView(index_input);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> b = 
    result_x.subView(index_dp);

  // compute (Jn)_x[A b]
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> tmp = 
    result_x.clone(NOX::ShapeCopy);
  status = group->computeDJnDxaMulti(*nullVector, *JnVector, result_x,
				     *tmp);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // compute (Jn)_x[A b] - [G d(Jn)/dp]
  tmp->update(-1.0, input_null, 1.0);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute [C d] = J^-1 (Jn)_x[A b] - [G d(Jn)/dp]
  status = group->applyJacobianInverseMultiVector(params, *tmp, result_null);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> C = 
    result_null.subView(index_input);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> d = 
    result_null.subView(index_dp);

  // compute z = (h + phi^T C) / phi^T d
  tpGroup->lTransNorm(*C, result_param);
  result_param += input_param;
  double denom = tpGroup->lTransNorm((*d)[0]);
  result_param.scale(1.0/denom);

  // compute A = A - b*z (remember A is a sub-view of result_x)
  A->update(Teuchos::NO_TRANS, -1.0, *b, result_param, 1.0);

  // compute C = -C + d*z (remember C is a sub-view of result_null)
  C->update(Teuchos::NO_TRANS, 1.0, *d, result_param, -1.0);

  return finalStatus;
}
