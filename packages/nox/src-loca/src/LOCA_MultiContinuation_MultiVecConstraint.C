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

#include "LOCA_MultiContinuation_MultiVecConstraint.H"

LOCA::MultiContinuation::MultiVecConstraint::MultiVecConstraint(
    const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& dx_) :
  dx(dx_->clone(NOX::DeepCopy)),
  x(dx->clone(1)),
  constraints(dx->numVectors(), 1),
  isValidConstraints(false)
{
}

LOCA::MultiContinuation::MultiVecConstraint::MultiVecConstraint(
		  const LOCA::MultiContinuation::MultiVecConstraint& source, 
		  NOX::CopyType type) : 
  dx(source.dx->clone(type)),
  x(source.x->clone(type)),
  constraints(source.constraints),
  isValidConstraints(false)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

LOCA::MultiContinuation::MultiVecConstraint::~MultiVecConstraint()
{
}

void
LOCA::MultiContinuation::MultiVecConstraint::setDx(
	    const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& dx_)
{
  *dx = *dx_;
}

void
LOCA::MultiContinuation::MultiVecConstraint::copy(
		   const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::MultiContinuation::MultiVecConstraint& source = 
    dynamic_cast<const LOCA::MultiContinuation::MultiVecConstraint&>(src);

  if (this != &source) {
    *dx = *source.dx;
    *x = *source.x;
    constraints.assign(source.constraints);
    isValidConstraints = source.isValidConstraints;
  }
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::MultiVecConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new MultiVecConstraint(*this, type));
}

int
LOCA::MultiContinuation::MultiVecConstraint::numConstraints() const
{
  return constraints.numRows();
}

void
LOCA::MultiContinuation::MultiVecConstraint::setX(
					      const NOX::Abstract::Vector& y)
{
  (*x)[0] = y;
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::MultiVecConstraint::setParam(int paramID, double val)
{
}

void
LOCA::MultiContinuation::MultiVecConstraint::setParams(
			 const vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::MultiVecConstraint::computeConstraints()
{
  if (!isValidConstraints) {
    x->multiply(1.0, *dx, constraints);
    isValidConstraints = true;
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::MultiVecConstraint::computeDX()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::MultiVecConstraint::computeDP(
		                const vector<int>& paramIDs, 
		                NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
				bool isValidG)
{
   string callingFunction = 
    "LOCA::MultiContinuation::MultiVecConstraint::computeDP()";
  NOX::Abstract::Group::ReturnType status;
  
  // Compute constraints if necessary
  if (!isValidG && !isValidConstraints)
    status = computeConstraints();

  if (!isValidG) {
    for (int i=0; i<constraints.numRows(); i++)
      dgdp(i,0) = constraints(i,0);
  }

  // Set rest of dgdp to zero
  for (unsigned int j=0; j<paramIDs.size(); j++)
    for (int i=0; i<constraints.numRows(); i++)
      dgdp(i,j+1) = 0.0;

  return status;
}

bool
LOCA::MultiContinuation::MultiVecConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::MultiVecConstraint::isDX() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::MultiVecConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::MultiVecConstraint::getDX() const
{
  return dx.get();
}

bool
LOCA::MultiContinuation::MultiVecConstraint::isDXZero() const
{
  return false;
}

void
LOCA::MultiContinuation::MultiVecConstraint::notifyCompletedStep()
{
}
