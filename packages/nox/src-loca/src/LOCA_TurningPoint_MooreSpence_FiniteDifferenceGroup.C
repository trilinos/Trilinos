// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"

LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup()
{
}

LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup(
         const LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup& source,
     NOX::CopyType type)
  : LOCA::MultiContinuation::FiniteDifferenceGroup(source, type)
{
}


LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::~FiniteDifferenceGroup()
{
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDpMulti(
                      const std::vector<int>& paramIDs,
                      const NOX::Abstract::Vector& nullVector,
                      NOX::Abstract::MultiVector& result,
                      bool isValid)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDp(*this, paramIDs, nullVector, result, isValid);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDxaMulti(
                   const NOX::Abstract::Vector& nullVector,
                   const NOX::Abstract::MultiVector& aVector,
                   NOX::Abstract::MultiVector& result)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDxa(*this, nullVector, aVector, result);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDxaMulti(
                   const NOX::Abstract::Vector& nullVector,
                   const NOX::Abstract::Vector& JnVector,
                   const NOX::Abstract::MultiVector& aVector,
                   NOX::Abstract::MultiVector& result)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDxa(*this, nullVector, aVector, JnVector, result);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDwtJnDxMulti(
                       const NOX::Abstract::MultiVector& w,
                       const NOX::Abstract::Vector& nullVector,
                       NOX::Abstract::MultiVector& result)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDwtJnDx(*this, w, nullVector, result);
}
