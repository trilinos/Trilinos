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

#include "LOCA_MultiContinuation_NaturalConstraint.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint(
    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
    const Teuchos::RefCountPtr<LOCA::MultiContinuation::NaturalGroup>& grp) :
  globalData(global_data),
  naturalGroup(grp),
  constraints(grp->getNumParams(), 1),
  dgdp(grp->getNumParams(), grp->getNumParams()),
  isValidConstraints(false)
{
  // Initialize dgdp to identity matrix
  dgdp.putScalar(0.0);
  for (int i=0; i<grp->getNumParams(); i++)
    dgdp(i,i) = 1.0;
}

LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint(
		  const LOCA::MultiContinuation::NaturalConstraint& source, 
		  NOX::CopyType type) : 
  globalData(source.globalData),
  naturalGroup(),
  constraints(source.constraints),
  dgdp(source.dgdp),
  isValidConstraints(source.isValidConstraints)
{
}

LOCA::MultiContinuation::NaturalConstraint::~NaturalConstraint()
{
}

LOCA::MultiContinuation::NaturalConstraint&
LOCA::MultiContinuation::NaturalConstraint::operator=(
		   const LOCA::MultiContinuation::NaturalConstraint& source)
{
  if (this != &source) {
    globalData = source.globalData;
    constraints.assign(source.constraints);
    dgdp.assign(source.dgdp);
    isValidConstraints = source.isValidConstraints;
  }

  return *this;
}

void
LOCA::MultiContinuation::NaturalConstraint::setNaturalGroup(const Teuchos::RefCountPtr<LOCA::MultiContinuation::NaturalGroup>& grp)
{
  naturalGroup = grp;
}

LOCA::MultiContinuation::ConstraintInterface& 
LOCA::MultiContinuation::NaturalConstraint::operator=(
		   const LOCA::MultiContinuation::ConstraintInterface& source)
{
  return operator=(dynamic_cast<const LOCA::MultiContinuation::NaturalConstraint&>(source));
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::NaturalConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new NaturalConstraint(*this, type));
}

int
LOCA::MultiContinuation::NaturalConstraint::numConstraints() const
{
  return constraints.numRows();
}

void
LOCA::MultiContinuation::NaturalConstraint::setX(
						const NOX::Abstract::Vector& y)
{
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::NaturalConstraint::setParam(int paramID, double val)
{
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::NaturalConstraint::setParams(
			 const vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  isValidConstraints = false;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  // Get current, previous solution vectors
  const LOCA::MultiContinuation::ExtendedVector& xVec = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(naturalGroup->
								 getX());
  const LOCA::MultiContinuation::ExtendedVector& prevXVec = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(naturalGroup->
								 getPrevX());

  for (int i=0; i<naturalGroup->getNumParams(); i++)
    constraints(i,0) = xVec.getScalar(i) - prevXVec.getScalar(i) - 
      naturalGroup->getStepSize(i);

  isValidConstraints = true;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::computeConstraintDerivatives()
{
  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiContinuation::NaturalConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::NaturalConstraint::isConstraintDerivatives() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::NaturalConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::NaturalConstraint::getConstraintDerivativesP() const
{
  return dgdp;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::applyConstraintDerivativesX(
		      double alpha, 
		      const NOX::Abstract::MultiVector& input_x,
		      NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{
  result_p.putScalar(0.0);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::applyConstraintDerivativesX(
		             Teuchos::ETransp transb,
			     double alpha, 
			     const NOX::Abstract::MultiVector::DenseMatrix& b,
			     double beta,
			     NOX::Abstract::MultiVector& result_x) const
{
  for (int i=0; i<result_x.numVectors(); i++)
    result_x[i].scale(beta);

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiContinuation::NaturalConstraint::isConstraintDerivativesXZero() const
{
  return true;
}

bool
LOCA::MultiContinuation::NaturalConstraint::isConstraintDerivativesPZero() const
{
  return false;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::computeDgDp(
				const vector<int>& paramIDs, 
		                NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
				bool isValidG)
{
  string callingFunction = 
    "LOCA::MultiContinuation::NaturalConstraint::computeDgDp()";
  globalData->locaErrorCheck->throwError(callingFunction,
					 "Natural Constraint does not support derivatives with respect to non-continuation parameters!");

  return NOX::Abstract::Group::NotDefined;
}
