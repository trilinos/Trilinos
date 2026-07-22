// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_StandardParameterEntryXMLConverters.hpp"

namespace ROL {


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


} //namespace ROL

