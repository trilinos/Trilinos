// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_ParameterEntryXMLConverter.hpp"
#include "ROL_ParameterEntryXMLConverterDB.hpp"

namespace ROL {

/**
 * \brief Thrown when a parameter entry tag is missing it's value attribute.
 */
class NoValueAttributeException : public std::logic_error{
public:
  /**
   * \brief Constructs a NoValueAttributeExecption.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoValueAttributeException(const std::string& what_arg):std::logic_error(what_arg){}
};

ParameterEntry
ParameterEntryXMLConverter::fromXMLtoParameterEntry(
  const XMLObject &xmlObj) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !xmlObj.hasAttribute(getValueAttributeName()),
    NoValueAttributeException,
    ParameterEntry::getTagName() <<" tags must "
    "have a " << getValueAttributeName() << " attribute" << std::endl <<
    "Bad Parameter: " <<
    xmlObj.getAttribute(getNameAttributeName()) <<
    std::endl << std::endl);

  ParameterEntry toReturn;
  bool isDefault = false;
  bool isUsed = false;
  std::string docString = "";


  if(xmlObj.hasAttribute(getDefaultAttributeName())){
    isDefault = xmlObj.getRequiredBool(getDefaultAttributeName());
  }

  if(xmlObj.hasAttribute(getUsedAttributeName())){
    isUsed = xmlObj.getRequiredBool(getUsedAttributeName());
  }

  if(xmlObj.hasAttribute(getDocStringAttributeName())){
    docString = xmlObj.getRequired(getDocStringAttributeName());
  }

  toReturn.setAnyValue(getAny(xmlObj), isDefault);
  toReturn.setDocString(docString);

  if(isUsed){
    toReturn.getAny();
  }

  return toReturn;
}

} // namespace ROL

