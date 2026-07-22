// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConditionXMLConverter.hpp"


namespace Teuchos{


RCP<Condition>
ConditionXMLConverter::fromXMLtoCondition(const XMLObject &xmlObj,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap) const
{
  return convertXML(xmlObj, entryIDsMap);
}

XMLObject
ConditionXMLConverter::fromConditiontoXML(
  const RCP<const Condition> condition,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const
{
  XMLObject toReturn(Condition::getXMLTagName());
  toReturn.addAttribute(
    getTypeAttributeName(), condition->getTypeAttributeValue());
  convertCondition(condition, toReturn, entryIDsMap);
  return toReturn;
}


} // namespace Teuchos

