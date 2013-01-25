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

#include "LOCA_MultiContinuation_NaturalConstraint.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::MultiContinuation::NaturalGroup>& grp) :
  globalData(global_data),
  naturalGroup(grp),
  constraints(grp->getNumParams(), 1),
  isValidConstraints(false),
  conParamIDs(grp->getContinuationParameterIDs())
{
}

LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint(
		  const LOCA::MultiContinuation::NaturalConstraint& source, 
		  NOX::CopyType type) : 
  globalData(source.globalData),
  naturalGroup(),
  constraints(source.constraints),
  isValidConstraints(source.isValidConstraints),
  conParamIDs(source.conParamIDs)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

LOCA::MultiContinuation::NaturalConstraint::~NaturalConstraint()
{
}

void
LOCA::MultiContinuation::NaturalConstraint::setNaturalGroup(const Teuchos::RCP<LOCA::MultiContinuation::NaturalGroup>& grp)
{
  naturalGroup = grp;
}

void
LOCA::MultiContinuation::NaturalConstraint::copy(
		   const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::MultiContinuation::NaturalConstraint& source = 
    dynamic_cast<const LOCA::MultiContinuation::NaturalConstraint&>(src);

  if (this != &source) {
    globalData = source.globalData;
    constraints.assign(source.constraints);
    isValidConstraints = source.isValidConstraints;
    conParamIDs = source.conParamIDs;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
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
			 const std::vector<int>& paramIDs, 
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
LOCA::MultiContinuation::NaturalConstraint::computeDX()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::computeDP(
		                const std::vector<int>& paramIDs, 
		                NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
				bool isValidG)
{
   std::string callingFunction = 
    "LOCA::MultiContinuation::NaturalConstraint::computeDP()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  
  // Compute constraints if necessary
  if (!isValidG && !isValidConstraints) {
    status = computeConstraints();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  if (!isValidG) {
    for (int i=0; i<constraints.numRows(); i++)
      dgdp(i,0) = constraints(i,0);
  }

  // If a param ID is equal to a constraint param ID, then that column
  // of dgdp is given by that column of the identity matrix, other wise
  // that column is zero
  std::vector<int>::const_iterator it;
  for (unsigned int i=0; i<paramIDs.size(); i++) {
    for (int k=0; k<constraints.numRows(); k++)
	dgdp(k,i+1) = 0.0;
    it = find(conParamIDs.begin(), conParamIDs.end(), paramIDs[i]);
    if (it != conParamIDs.end())
      dgdp(it-conParamIDs.begin(),i+1) = 1.0;
  }

  return finalStatus;
}

bool
LOCA::MultiContinuation::NaturalConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::NaturalConstraint::isDX() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::NaturalConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::NaturalConstraint::getDX() const
{
  return NULL;
}

bool
LOCA::MultiContinuation::NaturalConstraint::isDXZero() const
{
  return true;
}

