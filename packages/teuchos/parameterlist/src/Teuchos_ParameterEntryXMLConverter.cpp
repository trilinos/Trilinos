// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterEntryXMLConverter.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_ValidatorXMLConverter.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"

namespace Teuchos{


ParameterEntry
ParameterEntryXMLConverter::fromXMLtoParameterEntry(
  const XMLObject &xmlObj) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      xmlObj.getRequired(getTypeAttributeName()) != getTypeAttributeValue(),
      BadParameterEntryXMLConverterTypeException,
      "Error: this Parameter Entry XML tag has a type different than "
      "the XMLConverter being used to convert it." <<std::endl <<
      "Parameter name: " << xmlObj.getRequired(
      XMLParameterListWriter::getNameAttributeName()) << std::endl <<
      "XML Parameter Entry type: " <<
      xmlObj.getRequired(getTypeAttributeName()) << std::endl <<
      "XMLConverter type: " << getTypeAttributeValue() <<
      std::endl <<std::endl);
  #endif

  TEUCHOS_TEST_FOR_EXCEPTION(
    !xmlObj.hasAttribute(getValueAttributeName()),
    NoValueAttributeExecption,
    ParameterEntry::getTagName() <<" tags must "
    "have a " << getValueAttributeName() << " attribute" << std::endl <<
    "Bad Parameter: " <<
    xmlObj.getAttribute(XMLParameterListWriter::getNameAttributeName()) <<
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


XMLObject
ParameterEntryXMLConverter::fromParameterEntrytoXML(
  RCP<const ParameterEntry> entry,
  const std::string &name,
  const ParameterEntry::ParameterEntryID& id,
  const ValidatortoIDMap& validatorIDsMap) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    (entry->getAny().typeName() != getTypeAttributeValue())
    &&
    (
      getTypeAttributeValue() !=
      ParameterEntryXMLConverterDB::getDefaultConverter()->getTypeAttributeValue()
    ),
    BadParameterEntryXMLConverterTypeException,
    "Error: This converter can't convert the given ParameterEntry to XML "
    "because their types don't match." << std::endl <<
    "Parameter name: " << name << std::endl <<
    "Parameter type: " << entry->getAny().typeName() << std::endl <<
    "Converter type: " << getTypeAttributeValue() << std::endl << std::endl);
  #endif

  XMLObject toReturn(ParameterEntry::getTagName());
  toReturn.addAttribute(
    XMLParameterListWriter::getNameAttributeName(), name);
  toReturn.addAttribute(getTypeAttributeName(), getTypeAttributeValue());
  toReturn.addAttribute(getDocStringAttributeName(), entry->docString());
  toReturn.addAttribute(getIdAttributeName(), id);
  toReturn.addAttribute(
    getValueAttributeName(), getValueAttributeValue(entry));
  toReturn.addBool(getDefaultAttributeName(), entry->isDefault());
  toReturn.addBool(getUsedAttributeName(), entry->isUsed());
  if(nonnull(entry->validator())){
    TEUCHOS_TEST_FOR_EXCEPTION(
      validatorIDsMap.find(entry->validator()) == validatorIDsMap.end(),
      MissingValidatorDefinitionException,
      "Could not find validator in given ValidatorIDsMap! " <<
      std::endl << std::endl);
    toReturn.addAttribute(
      ValidatorXMLConverter::getIdAttributeName(),
      validatorIDsMap.find(entry->validator())->second);
  }
  return toReturn;
}


} // namespace Teuchos

