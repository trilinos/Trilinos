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

#include "LOCA_MultiContinuation_CompositeConstraint.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::MultiContinuation::CompositeConstraint::CompositeConstraint() :
  globalData(),
  numConstraintObjects(0),
  constraintPtrs(),
  indices(),
  totalNumConstraints(0),
  constraints(),
  isValidConstraints(false),
  isValidDX(false)
{
}

LOCA::MultiContinuation::CompositeConstraint::CompositeConstraint(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const std::vector< Teuchos::RCP<
    LOCA::MultiContinuation::ConstraintInterface> >& constraintObjects) :
  globalData(),
  numConstraintObjects(0),
  constraintPtrs(),
  indices(),
  totalNumConstraints(0),
  constraints(),
  isValidConstraints(false),
  isValidDX(false)
{
  init(global_data, constraintObjects);
}

LOCA::MultiContinuation::CompositeConstraint::CompositeConstraint(
		  const LOCA::MultiContinuation::CompositeConstraint& source, 
		  NOX::CopyType type) : 
  globalData(source.globalData),
  numConstraintObjects(source.numConstraintObjects),
  constraintPtrs(source.constraintPtrs),
  indices(source.indices),
  totalNumConstraints(source.totalNumConstraints),
  constraints(source.constraints),
  isValidConstraints(source.isValidConstraints),
  isValidDX(source.isValidDX)
{
}

LOCA::MultiContinuation::CompositeConstraint::~CompositeConstraint()
{
}

void
LOCA::MultiContinuation::CompositeConstraint::copy(
		   const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::MultiContinuation::CompositeConstraint& source = 
    dynamic_cast<const LOCA::MultiContinuation::CompositeConstraint&>(src);


  // I should copy the constraint objects instead of copying pointers?
  if (this != &source) {
    globalData = source.globalData;
    numConstraintObjects = source.numConstraintObjects;
    constraintPtrs = source.constraintPtrs;
    indices = source.indices;
    totalNumConstraints = source.totalNumConstraints;
    constraints.assign(source.constraints);
    isValidConstraints = source.isValidConstraints;
    isValidDX = source.isValidDX;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::CompositeConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new CompositeConstraint(*this, type));
}

int
LOCA::MultiContinuation::CompositeConstraint::numConstraints() const
{
  return totalNumConstraints;
}

void
LOCA::MultiContinuation::CompositeConstraint::setX(
					      const NOX::Abstract::Vector& y)
{
  for (int i=0; i<numConstraintObjects; i++)
    constraintPtrs[i]->setX(y);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::MultiContinuation::CompositeConstraint::setParam(int paramID, double val)
{
  for (int i=0; i<numConstraintObjects; i++)
    constraintPtrs[i]->setParam(paramID, val);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::MultiContinuation::CompositeConstraint::setParams(
			 const std::vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (int i=0; i<numConstraintObjects; i++)
    constraintPtrs[i]->setParams(paramIDs, vals);
  isValidConstraints = false;
  isValidDX = false;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::MultiContinuation::CompositeConstraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  const NOX::Abstract::MultiVector::DenseMatrix *g;

  for (int i=0; i<numConstraintObjects; i++) {
    status = constraintPtrs[i]->computeConstraints();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    g = &(constraintPtrs[i]->getConstraints());
    for (int j=0; j<constraintPtrs[i]->numConstraints(); j++)
      constraints(indices[i][j],0) = (*g)(j,0);
  }

  isValidConstraints = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::computeDX()
{
  if (isValidDX)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::MultiContinuation::CompositeConstraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<numConstraintObjects; i++) {
    status = constraintPtrs[i]->computeDX();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::computeDP(
		                const std::vector<int>& paramIDs, 
		                NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
				bool isValidG)
{
   std::string callingFunction = 
    "LOCA::MultiContinuation::CompositeConstraint::computeDP()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> dgdp_sub;
  int num_rows;
  int num_cols = dgdp.numCols();
  for (int i=0; i<numConstraintObjects; i++) {

    // Create a sub view of rows indices[i][0] -- indices[i][end] of dgdp
    num_rows = indices[i][constraintPtrs[i]->numConstraints()-1] - 
      indices[i][0] + 1;
    dgdp_sub = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
							       dgdp,
							       num_rows,
							       num_cols,
							       indices[i][0],
							       0));

    status = constraintPtrs[i]->computeDP(paramIDs, *dgdp_sub, isValidG);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}

bool
LOCA::MultiContinuation::CompositeConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::CompositeConstraint::isDX() const
{
  return isValidDX;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::CompositeConstraint::getConstraints() const
{
  return constraints;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::multiplyDX(
		      double alpha, 
		      const NOX::Abstract::MultiVector& input_x,
	              NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{
  std::string callingFunction = 
    "LOCA::MultiContinuation::CompositeConstraint::multiplyDX()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // If dg/dx is zero for every constraint, result_p is zero
  if (isDXZero()) {
    result_p.putScalar(0.0);
    return finalStatus;
  }

  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_p_sub;
  int num_rows;
  int num_cols = result_p.numCols();
  for (int i=0; i<numConstraintObjects; i++) {

    num_rows = constraintPtrs[i]->numConstraints();

    // if dg/dx is zero for this constraint, set corresponding entries of
    // result_p to zero
    if (constraintPtrs[i]->isDXZero()) {
      for (int j=0; j<num_rows; j++)
	for (int k=0; k<num_cols; k++)
	  result_p(indices[i][j],k) = 0.0;
    }
    else {

      // Create a sub view of rows indices[i][0] -- indices[i][end] 
      // of result_p 
      result_p_sub = 
	Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
								 result_p,
								 num_rows,
								 num_cols,
								 indices[i][0],
								 0));

      status = constraintPtrs[i]->multiplyDX(alpha, input_x, 
					     *result_p_sub);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
    }

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::addDX(
		              Teuchos::ETransp transb,
			      double alpha, 
	                      const NOX::Abstract::MultiVector::DenseMatrix& b,
			      double beta,
			      NOX::Abstract::MultiVector& result_x) const
{
  std::string callingFunction = 
    "LOCA::MultiContinuation::CompositeConstraint::addDX()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // First scale result_x
  result_x.scale(beta);

  // If dg/dx is zero for every constraint, result_x = beta * result_x
  if (isDXZero())
    return finalStatus;

  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> b_sub;
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x_sub;
  int num_rows;
  int num_cols = result_x.numVectors();
  for (int j=0; j<numConstraintObjects; j++) {
      
    if (!constraintPtrs[j]->isDXZero()) {

      // Create a sub view of rows indices[j][0] -- indices[j][end],
      // of b or b^T
      num_rows = constraintPtrs[j]->numConstraints();
      if (transb == Teuchos::NO_TRANS) {
	b_sub = 
	  Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							       Teuchos::View,
							       b,
							       num_rows,
							       num_cols,
							       indices[j][0],
							       0));
      }
      else {
	b_sub = 
	  Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(
							       Teuchos::View,
							       b,
							       num_cols,
							       num_rows,
							       0,
							       indices[j][0]));
      }
	
      // Compute result_x = result_x + alpha*(dg/dx)_j*op(b)_j
      status = constraintPtrs[j]->addDX(transb, alpha, *b_sub, 1.0, 
					result_x);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
    }

  }

  return finalStatus;
}

bool
LOCA::MultiContinuation::CompositeConstraint::isDXZero() const
{
  for (int i=0; i<numConstraintObjects; i++)
    if (!constraintPtrs[i]->isDXZero())
      return false;

  return true;
}

void
LOCA::MultiContinuation::CompositeConstraint::preProcessContinuationStep(
			   LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  for (int i=0; i<numConstraintObjects; i++)
    constraintPtrs[i]->preProcessContinuationStep(stepStatus);
}

void
LOCA::MultiContinuation::CompositeConstraint::postProcessContinuationStep(
			   LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  for (int i=0; i<numConstraintObjects; i++)
    constraintPtrs[i]->postProcessContinuationStep(stepStatus);
}

void
LOCA::MultiContinuation::CompositeConstraint::init(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const std::vector< Teuchos::RCP<
    LOCA::MultiContinuation::ConstraintInterface> >& constraintObjects)
{
  globalData = global_data;
  numConstraintObjects = constraintObjects.size();
  constraintPtrs = constraintObjects;
  indices.resize(numConstraintObjects);
  totalNumConstraints = 0;

  int n;
  for (int i=0; i<numConstraintObjects; i++) {
    n = constraintPtrs[i]->numConstraints();
    indices[i].resize(n);
    for (int j=0; j<n; j++)
      indices[i][j] = totalNumConstraints + j;
    totalNumConstraints += n;
  }
  constraints.shape(totalNumConstraints, 1);
}
