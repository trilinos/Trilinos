// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup()
{
}

LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup(
                const LOCA::MultiContinuation::FiniteDifferenceGroup& source,
        NOX::CopyType type)
{
  if (source.derivPtr != Teuchos::null)
    derivPtr = source.derivPtr->clone(type);
}

LOCA::MultiContinuation::FiniteDifferenceGroup::~FiniteDifferenceGroup()
{
}

void
LOCA::MultiContinuation::FiniteDifferenceGroup::copy(
                        const NOX::Abstract::Group& src)
{
  const LOCA::MultiContinuation::FiniteDifferenceGroup& source =
    dynamic_cast<const LOCA::MultiContinuation::FiniteDifferenceGroup&>(src);

  if (this != &source)
    if (source.derivPtr != Teuchos::null)
      derivPtr = source.derivPtr->clone();
}

NOX::Abstract::Group&
LOCA::MultiContinuation::FiniteDifferenceGroup::operator=(
                        const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

void
LOCA::MultiContinuation::FiniteDifferenceGroup::setDerivUtils(
              const Teuchos::RCP<LOCA::DerivUtils>& deriv)
{
  derivPtr = deriv;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::FiniteDifferenceGroup::computeDfDpMulti(
                      const std::vector<int>& paramIDs,
                      NOX::Abstract::MultiVector& dfdp,
                      bool isValidF)
{
  return derivPtr->computeDfDp(*this, paramIDs, dfdp, isValidF);
}

