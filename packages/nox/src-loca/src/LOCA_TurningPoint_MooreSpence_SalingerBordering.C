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

#include "LOCA_TurningPoint_MooreSpence_SalingerBordering.H"
#include "LOCA_TurningPoint_MooreSpence_ExtendedGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"

LOCA::TurningPoint::MooreSpence::SalingerBordering::SalingerBordering(
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
  dJndp()
{
}

LOCA::TurningPoint::MooreSpence::SalingerBordering::~SalingerBordering()
{
}

void
LOCA::TurningPoint::MooreSpence::SalingerBordering::setBlocks(
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
}

NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::SalingerBordering::solve(
	   Teuchos::ParameterList& params,
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& input,
           LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::SalingerBordering::solve()";
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
  
  // Create new multivectors with m+1 columns
  // First m columns store input_x, input_null, result_x, result_null
  // respectively, last column stores dfdp, dJndp, J^-1 dfdp, J^-1 dJndp
  // respectively
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+1);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_null = 
    input_null->clone(m+1);
  
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+1);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_null = 
    result_null->clone(m+1);
  
  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);

  // Set last column to dfdp
  (*cont_input_x)[m] = (*dfdp)[0];
  
    // Set first m columns to input_null
  cont_input_null->setBlock(*input_null, index_input);
  
  // Set last column to dJndp
  (*cont_input_null)[m] = (*dJndp)[0];
  
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
LOCA::TurningPoint::MooreSpence::SalingerBordering::solveTranspose(
	   Teuchos::ParameterList& params,
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& input,
           LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::SalingerBordering::solveTranspose()";
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
  
  // Create new multivectors with m+1 columns
  // First m columns store input_x, input_null, result_x, result_null
  // respectively, last column stores 0, -phi, J^-T tmp , -J^-T phi
  // respectively
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+1);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_null = 
    input_null->clone(m+1);
  
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+1);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_null = 
    result_null->clone(m+1);
  
  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);

  // Set last column to 0
  (*cont_input_x)[m].init(0);
  
    // Set first m columns to input_null
  cont_input_null->setBlock(*input_null, index_input);
  
  // Set last column to phi
  Teuchos::RCP<NOX::Abstract::Vector> phi = tpGroup->getLengthVector();
  (*cont_input_null)[m].update(-1.0, *phi, 0.0);
  
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
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::SalingerBordering::solveContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  int m = input_x.numVectors()-1;
  std::vector<int> index_input(m);
  std::vector<int> index_dp(1);
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
  Teuchos::RCP<NOX::Abstract::MultiVector> A = 
    result_x.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> b = 
    result_x.subView(index_dp);

  // compute (Jn)_x[A b]
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp = 
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
  Teuchos::RCP<NOX::Abstract::MultiVector> C = 
    result_null.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> d = 
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

// Solves turning point equations via classic Salinger bordering
// The first m columns of input_x and input_null store the RHS while
// the last column stores df/dp, d(Jn)/dp respectively.  Note however
// input_param has only m columns (not m+1).  result_x, result_null,
// are result_param have the same dimensions as their input counterparts
NOX::Abstract::Group::ReturnType 
LOCA::TurningPoint::MooreSpence::SalingerBordering::solveTransposeContiguous(
		  Teuchos::ParameterList& params,
		  const NOX::Abstract::MultiVector& input_x,
		  const NOX::Abstract::MultiVector& input_null,
	          const NOX::Abstract::MultiVector::DenseMatrix& input_param,
		  NOX::Abstract::MultiVector& result_x,
		  NOX::Abstract::MultiVector& result_null,
	          NOX::Abstract::MultiVector::DenseMatrix& result_param) const
{
  std::string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::SalingerBordering::solveTransposeContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast group to transpose solve group
  Teuchos::RCP<LOCA::Abstract::TransposeSolveGroup> tsGroup =
    Teuchos::rcp_dynamic_cast<LOCA::Abstract::TransposeSolveGroup>(group);
  if (tsGroup == Teuchos::null)
    globalData->locaErrorCheck->throwError(callingFunction,
					   "Underlying group must be derived from NOX::Abstract::TransposeSolveGroup for transpose solve");

  int m = input_x.numVectors()-1;
  std::vector<int> index_input(m);
  std::vector<int> index_dp(1);
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
  
  // compute [A b] = J^-T [G -phi]
  status = tsGroup->applyJacobianTransposeInverseMultiVector(params, 
							     input_null, 
							     result_null);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> A = 
    result_null.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> b = 
    result_null.subView(index_dp);

  // compute (Jn)_x^T[A b]
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp = 
    result_null.clone(NOX::ShapeCopy);
  status = group->computeDwtJnDxMulti(result_null, *nullVector, *tmp);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // compute [F 0] - (Jn)_x[A b]
  tmp->update(1.0, input_x, -1.0);

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute [C d] = J^-T( [F 0] - (Jn)_x[A b] )
  status = tsGroup->applyJacobianTransposeInverseMultiVector(params, *tmp, 
							     result_x);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> C = 
    result_x.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> d = 
    result_x.subView(index_dp);

  // compute (Jn)_p^T*[A b]
  NOX::Abstract::MultiVector::DenseMatrix t1(1,m+1);
  result_null.multiply(1.0, *dJndp, t1);

  // compute f_p^T*[C d]
  NOX::Abstract::MultiVector::DenseMatrix t2(1,m+1);
  result_x.multiply(1.0, *dfdp, t2);

  // compute z = (h - f_p^T C - (Jn)_p^TA) / (f_p^T*d + (Jn)_p^T*b)
  double denom = t2(0,m) + t1(0,m);
  for (int i=0; i<m; i++)
    result_param(0,i) = (input_param(0,i) - t2(0,i) - t1(0,i))/denom;

  // compute A = A + b*z (remember A is a sub-view of result_null)
  A->update(Teuchos::NO_TRANS, 1.0, *b, result_param, 1.0);

  // compute C = C + d*z (remember C is a sub-view of result_x)
  C->update(Teuchos::NO_TRANS, 1.0, *d, result_param, 1.0);

  return finalStatus;
}
