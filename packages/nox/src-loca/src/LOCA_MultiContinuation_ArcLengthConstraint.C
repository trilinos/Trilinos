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

#include "LOCA_MultiContinuation_ArcLengthConstraint.H"
#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint(
    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
    const Teuchos::RefCountPtr<LOCA::MultiContinuation::ArcLengthGroup>& grp) :
  globalData(global_data),
  arcLengthGroup(grp),
  constraints(grp->getNumParams(), 1),
  isValidConstraints(false)
{
}

LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint(
		  const LOCA::MultiContinuation::ArcLengthConstraint& source, 
		  NOX::CopyType type) : 
  globalData(source.globalData),
  arcLengthGroup(),
  constraints(source.constraints),
  isValidConstraints(source.isValidConstraints)
{
}

LOCA::MultiContinuation::ArcLengthConstraint::~ArcLengthConstraint()
{
}

LOCA::MultiContinuation::ArcLengthConstraint&
LOCA::MultiContinuation::ArcLengthConstraint::operator=(
		   const LOCA::MultiContinuation::ArcLengthConstraint& source)
{
  if (this != &source) {
    globalData = source.globalData;
    constraints.assign(source.constraints);
    isValidConstraints = source.isValidConstraints;
  }

  return *this;
}

void
LOCA::MultiContinuation::ArcLengthConstraint::setArcLengthGroup(const Teuchos::RefCountPtr<LOCA::MultiContinuation::ArcLengthGroup>& grp)
{
  arcLengthGroup = grp;
}

LOCA::MultiContinuation::ConstraintInterface& 
LOCA::MultiContinuation::ArcLengthConstraint::operator=(
		   const LOCA::MultiContinuation::ConstraintInterface& source)
{
  return operator=(dynamic_cast<const LOCA::MultiContinuation::ArcLengthConstraint&>(source));
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::ArcLengthConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new ArcLengthConstraint(*this, type));
}

int
LOCA::MultiContinuation::ArcLengthConstraint::numConstraints() const
{
  return constraints.numRows();
}

void
LOCA::MultiContinuation::ArcLengthConstraint::setX(
					      const NOX::Abstract::Vector& y)
{
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::ArcLengthConstraint::setParam(int paramID, double val)
{
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::ArcLengthConstraint::setParams(
			 const vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  isValidConstraints = false;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::MultiContinuation::ArcLengthConstraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute predictor if necessary
  if (!arcLengthGroup->isPredictor()) {
    status = arcLengthGroup->computePredictor();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Get tangent vector
  const LOCA::MultiContinuation::ExtendedMultiVector& scaledTangent = 
    arcLengthGroup->getScaledPredictorTangent();
  const LOCA::MultiContinuation::ExtendedMultiVector& tangent = 
    arcLengthGroup->getPredictorTangent();

  // Compute secant vector
  LOCA::MultiContinuation::ExtendedMultiVector* secant = 
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector*>(tangent.
								clone(1));
  (*secant)[0].update(1.0, arcLengthGroup->getX(), 
		      -1.0, arcLengthGroup->getPrevX(), 0.0);

  // Compute [dx/ds; dp/ds]^T * [x - x_o; p - p_o] - ds
  scaledTangent.multiply(1.0, *secant, constraints);
  for (int i=0; i<arcLengthGroup->getNumParams(); i++)
    constraints(i,0) -= arcLengthGroup->getStepSize(i) * 
      scaledTangent[i].dot(tangent[i]);

  delete secant;

  isValidConstraints = true;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeConstraintDerivatives()
{
  if (!isValidConstraints)
    return computeConstraints();
  else
    return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiContinuation::ArcLengthConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::ArcLengthConstraint::isConstraintDerivatives() const
{
  return isValidConstraints;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::ArcLengthConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector::DenseMatrix*
LOCA::MultiContinuation::ArcLengthConstraint::getConstraintDerivativesP() const
{
  const LOCA::MultiContinuation::ExtendedMultiVector& tangent = 
    arcLengthGroup->getScaledPredictorTangent();

  return &tangent.getScalars();
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::ArcLengthConstraint::getConstraintDerivativesX() const
{
  // Get tangent vector
  const LOCA::MultiContinuation::ExtendedMultiVector& tangent = 
    arcLengthGroup->getScaledPredictorTangent();

  return &tangent.getXMultiVec();
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::applyConstraintDerivativesX(
		      double alpha, 
		      const NOX::Abstract::MultiVector& input_x,
		      NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{
  // Get tangent vector
  const LOCA::MultiContinuation::ExtendedMultiVector& tangent = 
    arcLengthGroup->getScaledPredictorTangent();

  // Get x component of tangent
  const NOX::Abstract::MultiVector& tangent_x = tangent.getXMultiVec();

  // Multiply
  //tangent_x.multiply(alpha, input_x, result_p);
  input_x.multiply(alpha, tangent_x, result_p);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::applyConstraintDerivativesX(
		             Teuchos::ETransp transb,
			     double alpha, 
			     const NOX::Abstract::MultiVector::DenseMatrix& b,
			     double beta,
			     NOX::Abstract::MultiVector& result_x) const
{
  // Get tangent vector
  const LOCA::MultiContinuation::ExtendedMultiVector& tangent = 
    arcLengthGroup->getScaledPredictorTangent();

  // Get x component of tangent
  const NOX::Abstract::MultiVector& tangent_x = tangent.getXMultiVec();

  // Update
  result_x.update(transb, alpha, tangent_x, b, beta);

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiContinuation::ArcLengthConstraint::isConstraintDerivativesXZero() const
{
  return false;
}

bool
LOCA::MultiContinuation::ArcLengthConstraint::isConstraintDerivativesPZero() const
{
  return false;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeDgDp(
				const vector<int>& paramIDs, 
		                NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
				bool isValidG)
{
  string callingFunction = 
    "LOCA::MultiContinuation::ArcLengthConstraint::computeDgDp()";
  globalData->locaErrorCheck->throwError(callingFunction,
					 "ArcLength Constraint does not support derivatives with respect to non-continuation parameters!");

  return NOX::Abstract::Group::NotDefined;
}
