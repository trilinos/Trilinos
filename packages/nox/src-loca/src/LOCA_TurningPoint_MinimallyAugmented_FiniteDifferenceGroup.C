// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"

LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::
FiniteDifferenceGroup()
{
}

LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::
FiniteDifferenceGroup(
  const LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup& source,
  NOX::CopyType type)
  :  LOCA::MultiContinuation::FiniteDifferenceGroup(source, type),
     LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(source, type)
{
}


LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::
~FiniteDifferenceGroup()
{
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::
computeDwtJnDp(const std::vector<int>& paramIDs,
           const NOX::Abstract::Vector& w,
           const NOX::Abstract::Vector& nullVector,
           NOX::Abstract::MultiVector::DenseMatrix& result,
           bool isValid)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDwtJnDp(*this, paramIDs, w, nullVector, result, isValid);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::
computeDwtJDp(const std::vector<int>& paramIDs,
           const NOX::Abstract::Vector& w,
           NOX::Abstract::MultiVector& result,
           bool isValid)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDwtJDp(*this, paramIDs, w, result, isValid);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::
computeDwtJnDx(const NOX::Abstract::Vector& w,
           const NOX::Abstract::Vector& nullVector,
           NOX::Abstract::Vector& result)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDwtJnDx(*this, w, nullVector, result);
}
