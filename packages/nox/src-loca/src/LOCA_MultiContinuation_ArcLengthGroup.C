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

#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
				 LOCA::MultiContinuation::AbstractGroup& g, 
				 const vector<int>& paramIDs,
				 NOX::Parameter::List& params)
  : LOCA::MultiContinuation::ExtendedGroup(g, paramIDs, params),
    constraints(numParams, 1),
    isValidConstraints(false),
    isValidConstraintDerivatives(false),
    theta(numParams, 1.0)
{
}

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
			 const LOCA::MultiContinuation::ArcLengthGroup& source,
			 NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type),
    constraints(source.constraints),
    isValidConstraints(source.isValidConstraints),
    isValidConstraintDerivatives(source.isValidConstraintDerivatives),
    theta(source.theta)
{
}


LOCA::MultiContinuation::ArcLengthGroup::~ArcLengthGroup() 
{
}

LOCA::MultiContinuation::ArcLengthGroup&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
		       const LOCA::MultiContinuation::ArcLengthGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    LOCA::MultiContinuation::ExtendedGroup::operator=(source);
    constraints = source.constraints;
    isValidConstraints = source.isValidConstraints;
    isValidConstraintDerivatives = source.isValidConstraintDerivatives;
    theta = source.theta;
  }

  return *this;
}

LOCA::MultiContinuation::ExtendedGroup&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
		        const LOCA::MultiContinuation::ExtendedGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthGroup::computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::NultiContinuation::ArcLengthGroup::computeConstraintss()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // compute predictor
  if (!isValidConstraintDerivatives) {
    finalStatus = computeConstraintDerivatives();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }

  // create view of first column of xVec as a multivec
  LOCA::MultiContinuation::ExtendedMultiVector *xxMultiVec =
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector*>(xMultiVec.subView(index_f));

  // compute secant vector
  LOCA::MultiContinuation::ExtendedMultiVector *secantMultiVec =
    dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector*>(xxMultiVec->clone(NOX::DeepCopy));
  secantMultiVec->update(-1.0, prevXMultiVec, 1.0);
  
  computeScaledDotProduct(predictorMultiVec, *secantMultiVec, constraints);
  for (int i=0; i<numParams; i++)
    constraints(i,i) -= stepSize[i];

  delete secantMultiVec;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthGroup::computeConstraintDerivatives()
{
  if (isValidConstraintDerivatives)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::NultiContinuation::ArcLengthGroup::computeConstraintDerivatives()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute predictor
  finalStatus = computePredictor();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  isValidConstraintDerivatives = true;

  return finalStatus;
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraintDerivatives() const
{
  return isValidConstraintDerivatives;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::ArcLengthGroup::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::ArcLengthGroup::getConstraintDerivativesX() const
{
  return &(getPredictorDirections().getXMultiVec());
}

const NOX::Abstract::MultiVector::DenseMatrix*
LOCA::MultiContinuation::ArcLengthGroup::getConstraintDerivativesP() const
{
  return &(getPredictorDirections().getScalars());
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraintDerivativesXZero() const
{
  return false;
}

bool
LOCA::MultiContinuation::ArcLengthGroup::isConstraintDerivativesPZero() const
{
  return false;
}

LOCA::Extended::AbstractGroup&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
			      const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(source);
}

NOX::Abstract::Group* 
LOCA::MultiContinuation::ArcLengthGroup::clone(NOX::CopyType type) const
{
  return new LOCA::MultiContinuation::ArcLengthGroup(*this, type);
}

void
LOCA::MultiContinuation::ArcLengthGroup::setPrevX(
			    const LOCA::MultiContinuation::ExtendedVector& y) 
{
  LOCA::MultiContinuation::ExtendedGroup::setPrevX(y);
  isValidConstraints = false;
  isValidConstraintDerivatives = false;
}

void
LOCA::MultiContinuation::ArcLengthGroup::setStepSize(int i, double deltaS) 
{
  LOCA::MultiContinuation::ExtendedGroup::setStepSize(i, deltaS);
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::ArcLengthGroup::setScaleFactor(int i, double t) 
{
  theta[i] = t;
}

double
LOCA::MultiContinuation::ArcLengthGroup::getScaleFactor(int i) const {
  return theta[i];
}

void
LOCA::MultiContinuation::ArcLengthGroup::recalculateScaleFactor(
						   const vector<double>& dpds) 
{
}

void
LOCA::MultiContinuation::ArcLengthGroup::computeScaledDotProduct(
			 const NOX::Abstract::MultiVector& x,
			 const NOX::Abstract::MultiVector& y,
			 NOX::Abstract::MultiVector::DenseMatrix& result) const
{
  const LOCA::MultiContinuation::ExtendedMultiVector& mx = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(x);
  const LOCA::MultiContinuation::ExtendedMultiVector& my = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(y);

  my.multiply(1.0, mx, result);

  return;
}



void
LOCA::MultiContinuation::ArcLengthGroup::resetIsValid() {
  LOCA::MultiContinuation::ExtendedGroup::resetIsValid();
  isValidConstraints = false;
  isValidConstraintDerivatives = false;
}

