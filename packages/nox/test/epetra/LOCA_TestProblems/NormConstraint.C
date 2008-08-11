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

#include "NormConstraint.H"
#include "LOCA_Parameter_Vector.H"

NormConstraint::NormConstraint(int N, const LOCA::ParameterVector& pVec,
			       NOX::Abstract::Vector& cloneVec) :
  n(N),
  constraints(1,1),
  isValidConstraints(false),
  p(pVec),
  x(),
  isZeroDgDx(false)
{
  constraints.putScalar(0.0);
  x = cloneVec.createMultiVector(1);
}

NormConstraint::NormConstraint(const NormConstraint& source, 
			       NOX::CopyType type) :
  n(source.n),
  constraints(source.constraints),
  isValidConstraints(false),
  p(source.p),
  x(source.x->clone(type)),
  isZeroDgDx(source.isZeroDgDx)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

NormConstraint::~NormConstraint()
{
}

void
NormConstraint::copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const NormConstraint& source = dynamic_cast<const NormConstraint&>(src);

  if (this != &source) {
    n = source.n;
    constraints = source.constraints;
    isValidConstraints = source.isValidConstraints;
    p = source.p;
    *x = *source.x;
    isZeroDgDx = source.isZeroDgDx;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
NormConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new NormConstraint(*this, type));
}

int
NormConstraint::numConstraints() const
{
  return 1;
}

void
NormConstraint::setX(const NOX::Abstract::Vector& y)
{
  (*x)[0] = y;
  x->scale(1.0/n);
  isValidConstraints = false;
}

void
NormConstraint::setParam(int paramID, double val)
{
  p[paramID] = val;
  isValidConstraints = false;
}

void
NormConstraint::setParams(
			 const vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    p[paramIDs[i]] = vals(i,0);
  isValidConstraints = false;
}

NOX::Abstract::Group::ReturnType
NormConstraint::computeConstraints()
{
  constraints(0,0) = 0.5 * n * (*x)[0].innerProduct((*x)[0]) - p.getValue("Constraint Param");

  isValidConstraints = true;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
NormConstraint::computeDX()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
NormConstraint::computeDP(const vector<int>& paramIDs, 
			  NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
			  bool isValidG)
{
  if (!isValidG) {
    dgdp(0,0) = constraints(0,0);
  }

  for (unsigned int i=0; i<paramIDs.size(); i++) {
    if (p.getLabel(paramIDs[i]) == "Constraint Param") {
      dgdp(0,i+1) = -1.0;
    }
    else {
      dgdp(0,i+1) = 0.0;
    }
  }

  return NOX::Abstract::Group::Ok;
}

bool
NormConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
NormConstraint::isDX() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
NormConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
NormConstraint::getDX() const
{
  return x.get();
}

bool
NormConstraint::isDXZero() const
{
  return isZeroDgDx;
}

void
NormConstraint::setIsZeroDX(bool flag)
{
  isZeroDgDx = flag;
}

