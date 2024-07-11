// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"

LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::
FiniteDifferenceGroup()
{
}

LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::
FiniteDifferenceGroup(
  const LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup& source,
  NOX::CopyType type)
  :  LOCA::MultiContinuation::FiniteDifferenceGroup(source, type),
     LOCA::Hopf::MooreSpence::FiniteDifferenceGroup(source, type)
{
}


LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::
~FiniteDifferenceGroup()
{
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::
computeDwtCeDp(const std::vector<int>& paramIDs,
           const NOX::Abstract::Vector& w1,
           const NOX::Abstract::Vector& w2,
           const NOX::Abstract::Vector& y,
           const NOX::Abstract::Vector& z,
           double omega,
           NOX::Abstract::MultiVector::DenseMatrix& result_real,
           NOX::Abstract::MultiVector::DenseMatrix& result_imag,
           bool isValid)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDwtCeDp(*this, paramIDs, w1, w2, y, z, omega,
           result_real, result_imag, isValid);
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::
computeDwtCeDx(const NOX::Abstract::Vector& w1,
           const NOX::Abstract::Vector& w2,
           const NOX::Abstract::Vector& y,
           const NOX::Abstract::Vector& z,
           double omega,
           NOX::Abstract::Vector& result_real,
           NOX::Abstract::Vector& result_imag)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDwtCeDx(*this, w1, w2, y, z, omega, result_real, result_imag);
}
