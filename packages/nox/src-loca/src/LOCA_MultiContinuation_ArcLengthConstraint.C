// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_MultiContinuation_ArcLengthConstraint.H"
#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::MultiContinuation::ArcLengthGroup>& grp) :
  globalData(global_data),
  arcLengthGroup(grp),
  constraints(grp->getNumParams(), 1),
  isValidConstraints(false),
  conParamIDs(grp->getContinuationParameterIDs())
{
}

LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint(
          const LOCA::MultiContinuation::ArcLengthConstraint& source,
          NOX::CopyType type) :
  globalData(source.globalData),
  arcLengthGroup(),
  constraints(source.constraints),
  isValidConstraints(false),
  conParamIDs(source.conParamIDs)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

LOCA::MultiContinuation::ArcLengthConstraint::~ArcLengthConstraint()
{
}

void
LOCA::MultiContinuation::ArcLengthConstraint::setArcLengthGroup(const Teuchos::RCP<LOCA::MultiContinuation::ArcLengthGroup>& grp)
{
  arcLengthGroup = grp;
}

void
LOCA::MultiContinuation::ArcLengthConstraint::copy(
           const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::MultiContinuation::ArcLengthConstraint& source =
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthConstraint&>(src);

  if (this != &source) {
    globalData = source.globalData;
    constraints.assign(source.constraints);
    isValidConstraints = source.isValidConstraints;
    conParamIDs = source.conParamIDs;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
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
                          const NOX::Abstract::Vector& /* y */)
{
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::ArcLengthConstraint::setParam(int /* paramID */, double /* val */)
{
  isValidConstraints = false;
}

void
LOCA::MultiContinuation::ArcLengthConstraint::setParams(
             const std::vector<int>& /* paramIDs */,
             const NOX::Abstract::MultiVector::DenseMatrix& /* vals */)
{
  isValidConstraints = false;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::MultiContinuation::ArcLengthConstraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute predictor if necessary
  if (!arcLengthGroup->isPredictor()) {
    status = arcLengthGroup->computePredictor();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }

  // Get tangent vector
  const LOCA::MultiContinuation::ExtendedMultiVector& scaledTangent =
    arcLengthGroup->getScaledPredictorTangent();
  const LOCA::MultiContinuation::ExtendedMultiVector& tangent =
    arcLengthGroup->getPredictorTangent();

  // Compute secant vector
  Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> secant =
    Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(tangent.clone(1));
  (*secant)[0].update(1.0, arcLengthGroup->getX(),
              -1.0, arcLengthGroup->getPrevX(), 0.0);

  // Compute [dx/ds; dp/ds]^T * [x - x_o; p - p_o] - ds
  secant->multiply(1.0, scaledTangent, constraints);
  for (int i=0; i<arcLengthGroup->getNumParams(); i++)
    constraints(i,0) -= arcLengthGroup->getStepSize(i) *
      scaledTangent[i].innerProduct(tangent[i]);

  isValidConstraints = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeDX()
{
  if (!isValidConstraints)
    return computeConstraints();
  else
    return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeDP(
                        const std::vector<int>& paramIDs,
                        NOX::Abstract::MultiVector::DenseMatrix& dgdp,
                bool isValidG)
{
   std::string callingFunction =
    "LOCA::MultiContinuation::ArcLengthConstraint::computeDP()";
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

  // Get tangent vector
  const LOCA::MultiContinuation::ExtendedMultiVector& scaledTangent =
    arcLengthGroup->getScaledPredictorTangent();

  // If a param ID is equal to a constraint param ID, then that column
  // of dgdp is given by that column of the scaled predictor, other wise
  // that column is zero
  std::vector<int>::const_iterator it;
  int idx;
  for (unsigned int i=0; i<paramIDs.size(); i++) {
    it = find(conParamIDs.begin(), conParamIDs.end(), paramIDs[i]);
    if (it == conParamIDs.end())
      for (int k=0; k<constraints.numRows(); k++)
    dgdp(k,i+1) = 0.0;
    else {
      idx = it - conParamIDs.begin();
      for (int k=0; k<constraints.numRows(); k++)
    dgdp(k,i+1) = scaledTangent.getScalar(k,idx);
    }
  }

  return finalStatus;
}

bool
LOCA::MultiContinuation::ArcLengthConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::MultiContinuation::ArcLengthConstraint::isDX() const
{
  return isValidConstraints;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::ArcLengthConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::ArcLengthConstraint::getDX() const
{
  // Get tangent vector
  const LOCA::MultiContinuation::ExtendedMultiVector& tangent =
    arcLengthGroup->getScaledPredictorTangent();

  return tangent.getXMultiVec().get();
}

bool
LOCA::MultiContinuation::ArcLengthConstraint::isDXZero() const
{
  return false;
}
