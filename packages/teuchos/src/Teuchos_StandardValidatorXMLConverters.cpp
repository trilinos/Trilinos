// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_StandardValidatorXMLConverters.hpp"

/*! \file Teuchos_StandardValidatorXMLConverters.hpp
*/
namespace Teuchos {

RCP<ParameterEntryValidator> AnyNumberValidatorXMLConverter::convertXML(
 const XMLObject& xmlObj, 
 IDtoValidatorMap& validatorMap) const
{
  AnyNumberParameterEntryValidator dummyValidator;
  TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator.getXMLTagName(), 
    std::runtime_error, 
    "Cannot convert xmlObject to StringToIntegralValidator. Expected a " << dummyValidator.getXMLTagName() 
    << " tag but got a " << xmlObj.getTag() << "tag");
  AnyNumberParameterEntryValidator::AcceptedTypes acceptedTypes;
  acceptedTypes.allowInt(xmlObj.getRequiredBool(getAllowIntAttributeName()));
  acceptedTypes.allowDouble(xmlObj.getRequiredBool(getAllowDoubleAttributeName()));
  acceptedTypes.allowString(xmlObj.getRequiredBool(getAllowStringAttributeName()));
  return rcp(new AnyNumberParameterEntryValidator(AnyNumberParameterEntryValidator::getPrefferedTypeStringEnum(xmlObj.getRequired(getPrefferedTypeAttributeName())), acceptedTypes));
}

XMLObject AnyNumberValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator, 
  ValidatortoIDMap& validatorMap) const
{
  TEST_FOR_EXCEPTION(!isAppropriateConverter(validator), std::runtime_error, "An AnyNumberValidatorXMLConverter is not apporpriate for this type of validator.");
  RCP<const AnyNumberParameterEntryValidator> convertedValidator = rcp_static_cast<const AnyNumberParameterEntryValidator>(validator);
  XMLObject toReturn(validator->getXMLTagName());
  toReturn.addBool(getAllowIntAttributeName(), convertedValidator->isIntAllowed());
  toReturn.addBool(getAllowDoubleAttributeName(), convertedValidator->isDoubleAllowed());
  toReturn.addBool(getAllowStringAttributeName(), convertedValidator->isStringAllowed());
  toReturn.addAttribute(getPrefferedTypeAttributeName(), convertedValidator->getPrefferedTypeString(convertedValidator->getPreferredType()));
  return toReturn;
}

bool AnyNumberValidatorXMLConverter::isAppropriateConverter(
const RCP<const ParameterEntryValidator> validator) const
{
  return !(rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(validator).is_null());
}

RCP<ParameterEntryValidator> FileNameValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj, IDtoValidatorMap& validatorMap) const
{
  FileNameValidator dummyValidator;
  TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator.getXMLTagName(), 
    std::runtime_error, 
    "Cannot convert xmlObject to StringToIntegralValidator. Expected a " << dummyValidator.getXMLTagName() 
    << " tag but got a " << xmlObj.getTag() << "tag");
  return rcp(new FileNameValidator(xmlObj.getWithDefault<bool>(getFileMustExistAttributeName(), FileNameValidator::mustAlreadyExistDefault())));
}

XMLObject FileNameValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator, 
  ValidatortoIDMap& validatorMap) const
{
  TEST_FOR_EXCEPTION(!isAppropriateConverter(validator), std::runtime_error, "An FileNameValidatorXMLConverter is not apporpriate for this type of validator.");
  RCP<const FileNameValidator> convertedValidator = rcp_static_cast<const FileNameValidator>(validator);
  XMLObject toReturn(validator->getXMLTagName());
  toReturn.addBool(getFileMustExistAttributeName(), convertedValidator->fileMustExist());
  return toReturn;
}

bool FileNameValidatorXMLConverter::isAppropriateConverter(
  const RCP<const ParameterEntryValidator> validator) const
{
  return !(rcp_dynamic_cast<const FileNameValidator>(validator).is_null());
}

RCP<ParameterEntryValidator> StringValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj, IDtoValidatorMap& validatorMap) const
{
  StringValidator dummyValidator;
  TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator.getXMLTagName(), 
    std::runtime_error, 
    "Cannot convert xmlObject to StringValidator. Expected a " << dummyValidator.getXMLTagName() 
    << " tag but got a " << xmlObj.getTag() << "tag");
  if(xmlObj.numChildren()!=0){
    Array<std::string> strings(xmlObj.numChildren());
    for(int i=0; i<xmlObj.numChildren(); ++i){
      XMLObject currentChild = xmlObj.getChild(i);
      TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(), 
        std::runtime_error,  
        "Cannot convert xmlObject to StringToIntegralValidator." 
        << "\n Unrecognized tag: " << currentChild.getTag());
      strings[i] = (currentChild.getRequired(getStringValueAttributeName()));
    }
    return rcp(new StringValidator(strings));
  }
  return rcp(new StringValidator());
}

XMLObject StringValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  ValidatortoIDMap& validatorMap) const
{
  TEST_FOR_EXCEPTION(
    !isAppropriateConverter(validator), 
	std::runtime_error, 
	"An StringValidatorXMLConverter is not apporpriate for this type of validator.");
  XMLObject toReturn(validator->getXMLTagName());
  Array<std::string>::const_iterator it = validator->validStringValues()->begin();
  for(; it != validator->validStringValues()->end(); ++it){
    XMLObject stringTag(getStringTagName());
    stringTag.addAttribute(getStringValueAttributeName(), *it);
    toReturn.addChild(stringTag);
  }
  return toReturn;
}

bool StringValidatorXMLConverter::isAppropriateConverter(
  const RCP<const ParameterEntryValidator> validator) const
{
  return !(rcp_dynamic_cast<const StringValidator>(validator).is_null());
}

RCP<ParameterEntryValidator> UnknownValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj,
  IDtoValidatorMap& validatorMap) const
{
  throw std::runtime_error("Unknown xml tag. Can't convert to a Validator.");
  return Teuchos::null;
}

XMLObject UnknownValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator, 
  ValidatortoIDMap& validatorMap) const
{
  throw std::runtime_error("Unknown validator. Convert to XML.");
  return NULL;
}

bool UnknownValidatorXMLConverter::isAppropriateConverter(
  const RCP<const ParameterEntryValidator> validator) const
{
  return true;
}

}

