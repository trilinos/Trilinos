// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ValidatorXMLConverter.hpp"

namespace Teuchos{

RCP<ParameterEntryValidator>
ValidatorXMLConverter::fromXMLtoValidator(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& validatorIDsMap) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
  RCP<const ParameterEntryValidator> dummyValidator = getDummyValidator();
  TEUCHOS_TEST_FOR_EXCEPTION(
    xmlObj.getRequired(getTypeAttributeName())
    !=
    dummyValidator->getXMLTypeName(),
    BadValidatorXMLConverterException,
    "Cannot convert xmlObject " <<
    ". Expected a " << dummyValidator->getXMLTypeName() <<
    " tag but got a " << xmlObj.getRequired(getTypeAttributeName()) << "type");
  #endif
  RCP<ParameterEntryValidator> toReturn =
    convertXML(xmlObj, validatorIDsMap);
  return toReturn;
}

XMLObject
ValidatorXMLConverter::fromValidatortoXML(
  const RCP<const ParameterEntryValidator> validator,
  const ValidatortoIDMap& validatorIDsMap,
  bool assignedID) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
  RCP<const ParameterEntryValidator> dummyValidator = getDummyValidator();
  TEUCHOS_TEST_FOR_EXCEPTION(
    validator->getXMLTypeName()
    !=
    dummyValidator->getXMLTypeName(),
    BadValidatorXMLConverterException,
    "Cannot convert Validator " <<
    ". Expected a " << dummyValidator->getXMLTypeName() <<
    " validator but got a " << validator->getXMLTypeName() << "type");
  #endif
  XMLObject toReturn(getValidatorTagName());
  toReturn.addAttribute(getTypeAttributeName(), validator->getXMLTypeName());
  if(assignedID){
    TEUCHOS_TEST_FOR_EXCEPTION(validatorIDsMap.find(validator) == validatorIDsMap.end(),
      MissingValidatorDefinitionException,
      "Could not find an id associated with the validator in the "
      "given validatorIDsMap to use when " <<
      "writing it to XML!" << std::endl << std::endl);
    toReturn.addAttribute(getIdAttributeName(),
      validatorIDsMap.find(validator)->second);
  }
  convertValidator(validator, toReturn, validatorIDsMap);
  return toReturn;
}


}

