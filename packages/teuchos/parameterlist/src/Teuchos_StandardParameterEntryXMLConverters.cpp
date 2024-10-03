// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_StandardParameterEntryXMLConverters.hpp"

namespace Teuchos{


const std::string AnyParameterEntryConverter::getTypeAttributeValue() const{
  return "any";
}

const std::string AnyParameterEntryConverter::getValueAttributeValue(
  RCP<const ParameterEntry> entry) const
{
  return toString(entry->getAny(false));
}

any AnyParameterEntryConverter::getAny(const XMLObject& xmlObj) const {
  return any(xmlObj.getRequired(getValueAttributeName()));
}


} //namespace Teuchos

