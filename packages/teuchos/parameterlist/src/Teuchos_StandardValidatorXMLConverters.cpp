// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_StandardValidatorXMLConverters.hpp"
#include "Teuchos_RCP.hpp"

/*! \file Teuchos_StandardValidatorXMLConverters.hpp
*/


namespace Teuchos {

RCP<ParameterEntryValidator> BoolValidatorXMLConverter::convertXML(
  const XMLObject& /* xmlObj */,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  return boolParameterEntryValidator();
}


void BoolValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> /* validator */,
  XMLObject& /* xmlObj */,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  //RCP<const AnyNumberParameterEntryValidator> castedValidator =
  //  rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(validator, true);

  // currently no action
}

#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator>
BoolValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<BoolParameterEntryValidator>::getDummyObject();
}
#endif

RCP<ParameterEntryValidator> AnyNumberValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
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
      acceptedTypes);
}


void AnyNumberValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const AnyNumberParameterEntryValidator> castedValidator =
    rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(validator, true);
  xmlObj.addBool(
    getAllowIntAttributeName(), castedValidator->isIntAllowed());
  xmlObj.addBool(
    getAllowDoubleAttributeName(), castedValidator->isDoubleAllowed());
  xmlObj.addBool(
    getAllowStringAttributeName(), castedValidator->isStringAllowed());
  xmlObj.addAttribute(getPrefferedTypeAttributeName(),
    castedValidator->getPrefferedTypeString(
      castedValidator->getPreferredType()));
}

#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator>
AnyNumberValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<AnyNumberParameterEntryValidator>::getDummyObject();
}
#endif

RCP<ParameterEntryValidator> FileNameValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  return rcp(
    new FileNameValidator(
      xmlObj.getWithDefault<bool>(
        getFileMustExistAttributeName(),
        FileNameValidator::mustAlreadyExistDefault()
        )
      )
    );
}


void FileNameValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const FileNameValidator> castedValidator =
    rcp_dynamic_cast<const FileNameValidator>(validator);
  xmlObj.addBool(
    getFileMustExistAttributeName(), castedValidator->fileMustExist());
}


#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator>
FileNameValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<FileNameValidator>::getDummyObject();
}
#endif


RCP<ParameterEntryValidator> StringValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  Array<std::string> strings(xmlObj.numChildren());
  if(xmlObj.numChildren()!=0){
    for(int i=0; i<xmlObj.numChildren(); ++i){
      XMLObject currentChild = xmlObj.getChild(i);
      TEUCHOS_TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(),
        BadTagException,
        "Error converting xmlObject to StringValidator." << std::endl <<
		    "Unrecognized tag: " << currentChild.getTag());
        strings[i] = (currentChild.getRequired(getStringValueAttributeName()));
    }
  }
  return rcp(new StringValidator(strings));
}


void StringValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const StringValidator> castedValidator =
    rcp_dynamic_cast<const StringValidator>(validator);

  if(!is_null(validator->validStringValues())){
    Array<std::string>::const_iterator it =
     validator->validStringValues()->begin();
    for(; it != validator->validStringValues()->end(); ++it){
      XMLObject stringTag(getStringTagName());
      stringTag.addAttribute(getStringValueAttributeName(), *it);
      xmlObj.addChild(stringTag);
    }
  }
}


#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator>
StringValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<StringValidator>::getDummyObject();
}
#endif

} // namespace Teuchos

