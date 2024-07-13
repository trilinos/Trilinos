// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Extended_MultiAbstractGroup.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup() const
{
  // First get the underlying group
  Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup> ulg =
    getUnderlyingGroup();

  // Cast underlying group to an extended group
  Teuchos::RCP<const LOCA::Extended::MultiAbstractGroup> ulgPtr =
    Teuchos::rcp_dynamic_cast<const LOCA::Extended::MultiAbstractGroup>(ulg);

  if (ulgPtr.get() == NULL) {
    // Underlying group is not extended, therefore return it
    return ulg;
  }

  else {
    // Underlying group is extended, therefore return its baselevel group
    return ulgPtr->getBaseLevelUnderlyingGroup();
  }

}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup()
{
  // First get the underlying group
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> ulg =
    getUnderlyingGroup();

  // Cast underlying group to an extended group
  Teuchos::RCP<LOCA::Extended::MultiAbstractGroup> ulgPtr =
    Teuchos::rcp_dynamic_cast<LOCA::Extended::MultiAbstractGroup>(ulg);

  if (ulgPtr.get() == NULL) {
    // Underlying group is not extended, therefore return it
    return ulg;
  }

  else {
    // Underlying group is extended, therefore return its baselevel group
    return ulgPtr->getBaseLevelUnderlyingGroup();
  }

}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Extended::MultiAbstractGroup::getNestedGroup()
{
  return this->getUnderlyingGroup();
}

Teuchos::RCP<const NOX::Abstract::Group>
LOCA::Extended::MultiAbstractGroup::getNestedGroup() const
{
  return this->getUnderlyingGroup();
}
