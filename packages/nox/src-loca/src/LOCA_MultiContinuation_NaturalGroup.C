// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
				 LOCA::MultiContinuation::AbstractGroup& g, 
				 const vector<int>& paramIDs,
				 NOX::Parameter::List& params)
  : LOCA::MultiContinuation::ExtendedGroup(g, paramIDs, params),
    constraints(numParams, 1),
    dgdp(numParams, numParams),
    isValidConstraints(false)
{
  // Initialize dgdp to identity matrix
  dgdp.putScalar(0.0);
  for (int i=0; i<numParams; i++)
    dgdp(i,i) = 1.0;
}

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
				 LOCA::MultiContinuation::AbstractGroup& g, 
				 const string& paramID,
				 NOX::Parameter::List& params)
  : LOCA::MultiContinuation::ExtendedGroup(g, paramID, params),
    constraints(numParams, 1),
    dgdp(numParams, numParams),
    isValidConstraints(false)
{
  // Initialize dgdp to identity matrix
  dgdp.putScalar(0.0);
  for (int i=0; i<numParams; i++)
    dgdp(i,i) = 1.0;
}

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
			 const LOCA::MultiContinuation::NaturalGroup& source,
			 NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type),
    constraints(source.constraints),
    dgdp(source.dgdp),
    isValidConstraints(source.isValidConstraints)
{
}


LOCA::MultiContinuation::NaturalGroup::~NaturalGroup() 
{
}

LOCA::MultiContinuation::NaturalGroup&
LOCA::MultiContinuation::NaturalGroup::operator=(
		       const LOCA::MultiContinuation::NaturalGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    LOCA::MultiContinuation::ExtendedGroup::operator=(source);
    constraints.assign(source.constraints);
    dgdp.assign(source.dgdp);
    isValidConstraints = source.isValidConstraints;
  }

  return *this;
}

LOCA::MultiContinuation::ExtendedGroup&
LOCA::MultiContinuation::NaturalGroup::operator=(
		        const LOCA::MultiContinuation::ExtendedGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::NaturalGroup&>(source);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalGroup::computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  for (int i=0; i<numParams; i++) {
    constraints(i,0) = 
      xMultiVec.getScalar(i,0) - prevXMultiVec.getScalar(i,0) - stepSize[i];
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalGroup::computeConstraintDerivatives()
{
  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiContinuation::NaturalGroup::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::NaturalGroup::isConstraintDerivatives() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::NaturalGroup::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::NaturalGroup::getConstraintDerivativesX() const
{
  return NULL;
}

const NOX::Abstract::MultiVector::DenseMatrix*
LOCA::MultiContinuation::NaturalGroup::getConstraintDerivativesP() const
{
  return &(dgdp);
}

bool
LOCA::MultiContinuation::NaturalGroup::isConstraintDerivativesXZero() const
{
  return true;
}

bool
LOCA::MultiContinuation::NaturalGroup::isConstraintDerivativesPZero() const
{
  return false;
}

LOCA::Extended::AbstractGroup&
LOCA::MultiContinuation::NaturalGroup::operator=(
			      const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::NaturalGroup&>(source);
}

NOX::Abstract::Group&
LOCA::MultiContinuation::NaturalGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::NaturalGroup&>(source);
}

NOX::Abstract::Group* 
LOCA::MultiContinuation::NaturalGroup::clone(NOX::CopyType type) const
{
  return new LOCA::MultiContinuation::NaturalGroup(*this, type);
}

void
LOCA::MultiContinuation::NaturalGroup::setPrevX(
					      const NOX::Abstract::Vector& y) 
{
  LOCA::MultiContinuation::ExtendedGroup::setPrevX(y);
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::NaturalGroup::setStepSize(double deltaS, int i) 
{
  LOCA::MultiContinuation::ExtendedGroup::setStepSize(deltaS, i);
  isValidConstraints = false;
}


NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalGroup::applyJacobianInverseNewton(
						NOX::Parameter::List& params) 
{
  // This method is specialized to the Newton solve where the right-hand-side
  // is f.  We take advantage of the fact that f and df/dp are in a 
  // contiguous multivector

  string callingFunction = 
    "LOCA::MultiContinuation::NaturalGroup::applyJacobianInverseNewton()";
  
  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Get x, param components of f vector (we only want the parameter 
  // components of f, not df/dp)
  const NOX::Abstract::Vector& f_x = fVec->getXVec();
  const NOX::Abstract::MultiVector::DenseMatrix& f_p = fVec->getScalars();

  // Get references to x, param components of newton vector
  NOX::Abstract::Vector& newton_x = newtonVec->getXVec();
  NOX::Abstract::MultiVector::DenseMatrix& newton_p = newtonVec->getScalars();

  // Call bordered solver applyInverse method
  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianInverse(params, f_x, newton_x);
  newton_p.putScalar(0.0);

  return status;
}

void
LOCA::MultiContinuation::NaturalGroup::resetIsValid() {
  LOCA::MultiContinuation::ExtendedGroup::resetIsValid();
  isValidConstraints = false;
}

