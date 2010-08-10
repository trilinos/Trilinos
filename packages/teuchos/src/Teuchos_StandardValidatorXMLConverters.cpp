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
    ParameterEntryValidator::ValidatorID validatorID) const;
{

  AnyNumberParameterEntryValidator::AcceptedTypes acceptedTypes;
  acceptedTypes.allowInt(xmlObj.getRequiredBool(getAllowIntAttributeName()));
  acceptedTypes.allowDouble(
    xmlObj.getRequiredBool(getAllowDoubleAttributeName()));
  acceptedTypes.allowString(
    xmlObj.getRequiredBool(getAllowStringAttributeName()));
  return anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::getPrefferedTypeStringEnum(
        xmlObj.getRequired(getPrefferedTypeAttributeName())), 
      acceptedTypes,
      validatorID);
}


XMLObject AnyNumberValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator) const
{
  RCP<const AnyNumberParameterEntryValidator> castedValidator = 
    rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(validator, true);
  XMLObject toReturn(validator->getXMLTagName());
  toReturn.addBool(
    getAllowIntAttributeName(), castedValidator->isIntAllowed());
  toReturn.addBool(
    getAllowDoubleAttributeName(), castedValidator->isDoubleAllowed());
  toReturn.addBool(
    getAllowStringAttributeName(), castedValidator->isStringAllowed());
  toReturn.addAttribute(getPrefferedTypeAttributeName(),
    castedValidator->getPrefferedTypeString(
      castedValidator->getPreferredType()));
  return toReturn;
}

#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator> 
AnyNumberValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<AnyNumberParameterEntryValidator>::getDummyObject();
}
#endif

RCP<ParameterEntryValidator> FileNameValidatorXMLConverter::convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;
{
  return rcp(
    new FileNameValidator(
      xmlObj.getWithDefault<bool>(
        getFileMustExistAttributeName(),
        FileNameValidator::mustAlreadyExistDefault()
        ),
      validatorID
      )
    );
}


XMLObject FileNameValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator) const
{
  RCP<const FileNameValidator> castedValidator =
    rcp_dynamic_cast<const FileNameValidator>(validator);
  XMLObject toReturn(validator->getXMLTagName());
  toReturn.addBool(
    getFileMustExistAttributeName(), castedValidator->fileMustExist());
  return toReturn;
}


#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator> 
FileNameValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<FileNameValidator>::getDummyObject();
}
#endif


RCP<ParameterEntryValidator> StringValidatorXMLConverter::convertXML(
    const XMLObject& xmlObj,
    ParameterEntryValidator::ValidatorID validatorID) const;
{
  Array<std::string> strings(xmlObj.numChildren());
  if(xmlObj.numChildren()!=0){
    for(int i=0; i<xmlObj.numChildren(); ++i){
      XMLObject currentChild = xmlObj.getChild(i);
      TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(), 
        BadTagException,  
        "Error converting xmlObject to StringValidator." << std::endl << 
		    "Unrecognized tag: " << currentChild.getTag());
        strings[i] = (currentChild.getRequired(getStringValueAttributeName()));
    }
  }
  return rcp(new StringValidator(strings, ValidatorID));
}


XMLObject StringValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator) const
{
  RCP<const StringValidator> castedValidator = 
    rcp_dynamic_cast<const StringValidator>(validator);
  XMLObject toReturn(validator->getXMLTagName());
  Array<std::string>::const_iterator it = 
    validator->validStringValues()->begin();
  for(; it != validator->validStringValues()->end(); ++it){
    XMLObject stringTag(getStringTagName());
    stringTag.addAttribute(getStringValueAttributeName(), *it);
    toReturn.addChild(stringTag);
  }
  return toReturn;
}


#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator> 
StringValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<StringValidator>::getDummyObject();
}
#endif

} // namespace Teuchos

