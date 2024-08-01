// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"

LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup()
{
}

LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup(
            const LOCA::Hopf::MooreSpence::FiniteDifferenceGroup& source,
        NOX::CopyType type)
  : LOCA::MultiContinuation::FiniteDifferenceGroup(source, type),
    LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(source, type)
{
}


LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::~FiniteDifferenceGroup()
{
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDp(
                    const std::vector<int>& paramIDs,
                const NOX::Abstract::Vector& yVector,
                const NOX::Abstract::Vector& zVector,
                double w,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag,
                bool isValid)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDp(*this, paramIDs, yVector, zVector, w,
         result_real, result_imag, isValid);
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa(
                const NOX::Abstract::Vector& yVector,
                const NOX::Abstract::Vector& zVector,
                double w,
                const NOX::Abstract::MultiVector& aVector,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDxa(*this, yVector, zVector, w, aVector,
          result_real, result_imag);
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa(
                const NOX::Abstract::Vector& yVector,
                const NOX::Abstract::Vector& zVector,
                double w,
                const NOX::Abstract::MultiVector& aVector,
                const NOX::Abstract::Vector& Ce_real,
                const NOX::Abstract::Vector& Ce_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDxa(*this, yVector, zVector, w, aVector, Ce_real, Ce_imag,
          result_real, result_imag);
}
