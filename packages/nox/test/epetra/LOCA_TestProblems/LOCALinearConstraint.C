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

#include "LOCALinearConstraint.H"
#include "LOCA_Parameter_Vector.H"

LinearConstraint::LinearConstraint(int NumConstraints, 
				   const LOCA::ParameterVector& pVec,
				   NOX::Abstract::Vector& cloneVec) :
  m(NumConstraints),
  constraints(m,1),
  isValidConstraints(false),
  x(),
  dgdx(),
  isZeroDgDx(false),
  p(pVec),
  pvec(Teuchos::View, p.getDoubleArrayPointer(), p.length(), p.length(), 1),
  dgdp(NumConstraints, pVec.length())
{
  constraints.putScalar(0.0);
  x = cloneVec.createMultiVector(1);
  dgdx = cloneVec.createMultiVector(m);
  dgdp.putScalar(0.0);
}

LinearConstraint::LinearConstraint(const LinearConstraint& source, 
				   NOX::CopyType type) :
  m(source.m),
  constraints(source.constraints),
  isValidConstraints(false),
  x(source.x->clone(type)),
  dgdx(source.dgdx->clone(type)),
  isZeroDgDx(source.isZeroDgDx),
  p(source.p),
  pvec(source.pvec),
  dgdp(source.dgdp)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

LinearConstraint::~LinearConstraint()
{
}

void
LinearConstraint::copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LinearConstraint& source = dynamic_cast<const LinearConstraint&>(src);

  if (this != &source) {
    m = source.m;
    constraints = source.constraints;
    isValidConstraints = source.isValidConstraints;
    *x = *source.x;
    *dgdx = *source.dgdx;
    isZeroDgDx = source.isZeroDgDx;
    p = source.p;
    pvec.assign(source.pvec);
    dgdp = source.dgdp;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LinearConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LinearConstraint(*this, type));
}

int
LinearConstraint::numConstraints() const
{
  return m;
}

void
LinearConstraint::setX(const NOX::Abstract::Vector& y)
{
  (*x)[0] = y;
  isValidConstraints = false;
}

void
LinearConstraint::setParam(int paramID, double val)
{
  p[paramID] = val;
}

void
LinearConstraint::setParams(
			 const std::vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    p[paramIDs[i]] = vals(i,0);
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeConstraints()
{
  // compute dg/dx*x + dg/dp*p
  x->multiply(1.0, *dgdx, constraints);
  if (pvec.numRows() > 0)
    constraints.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, dgdp, 
			 pvec, 1.0);
  isValidConstraints = true;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeDX()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeDP(const std::vector<int>& paramIDs, 
			  NOX::Abstract::MultiVector::DenseMatrix& dp, 
			  bool isValidG)
{
  if (!isValidG) {
    for (int i=0; i<m; i++)
      dp(i,0) = constraints(i,0);
  }

  for (unsigned int i=0; i<paramIDs.size(); i++) {
    for (int j=0; j<m; j++)
      dp(j,i+1) = dgdp(j,paramIDs[i]);
  }

  return NOX::Abstract::Group::Ok;
}

bool
LinearConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LinearConstraint::isDX() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LinearConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LinearConstraint::getDX() const
{
  return dgdx.get();
}

bool
LinearConstraint::isDXZero() const
{
  return isZeroDgDx;
}

void
LinearConstraint::setDgDx(const NOX::Abstract::MultiVector& A)
{
  *dgdx = A;
  isValidConstraints = false;
}

void
LinearConstraint::setDgDp(const NOX::Abstract::MultiVector::DenseMatrix& A)
{
  dgdp.assign(A);
  isValidConstraints = false;
}

void
LinearConstraint::setIsZeroDX(bool flag)
{
  isZeroDgDx = flag;
}
