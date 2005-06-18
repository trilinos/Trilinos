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

#include "LinearConstraint.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Utils.H"

LinearConstraint::LinearConstraint(int NumConstraints, 
				   NOX::Abstract::Vector& cloneVec) :
  m(NumConstraints),
  constraints(m,1),
  isValidConstraints(false),
  x(),
  dgdx(),
  isZeroDgDx(false)
{
  constraints.putScalar(0.0);
  x = Teuchos::rcp(cloneVec.createMultiVector(1));
  dgdx = Teuchos::rcp(cloneVec.createMultiVector(m));
}

LinearConstraint::LinearConstraint(const LinearConstraint& source, 
				   NOX::CopyType type) :
  m(source.m),
  constraints(source.constraints),
  isValidConstraints(false),
  x(Teuchos::rcp(source.x->clone(type))),
  dgdx(Teuchos::rcp(source.dgdx->clone(type))),
  isZeroDgDx(source.isZeroDgDx)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

LinearConstraint::~LinearConstraint()
{
}

LinearConstraint&
LinearConstraint::operator=(const LinearConstraint& source)
{
  if (this != &source) {
    m = source.m;
    constraints = source.constraints;
    isValidConstraints = source.isValidConstraints;
    *x = *source.x;
    *dgdx = *source.dgdx;
    isZeroDgDx = source.isZeroDgDx;
  }

  return *this;
}

LOCA::MultiContinuation::ConstraintInterface& 
LinearConstraint::operator=(
		   const LOCA::MultiContinuation::ConstraintInterface& source)
{
  return operator=(dynamic_cast<const LinearConstraint&>(source));
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
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
}

void
LinearConstraint::setParams(
			 const vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeConstraints()
{
  x->multiply(1.0, *dgdx, constraints);
  isValidConstraints = true;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeDX()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeDP(const vector<int>& paramIDs, 
			  NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
			  bool isValidG)
{
  if (!isValidG) {
    for (int i=0; i<m; i++)
      dgdp(i,0) = constraints(i,0);
  }

  for (unsigned int i=0; i<paramIDs.size(); i++) {
    for (int j=0; j<m; j++)
      dgdp(j,i+1) = 0.0;
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
LinearConstraint::setIsZeroDX(bool flag)
{
  isZeroDgDx = flag;
}
