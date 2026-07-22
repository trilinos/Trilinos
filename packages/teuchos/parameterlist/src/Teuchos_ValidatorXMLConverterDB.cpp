// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_StaticSetupMacro.hpp"



namespace Teuchos {


void ValidatorXMLConverterDB::addConverter(
  RCP<const ParameterEntryValidator> validator,
  RCP<ValidatorXMLConverter> converterToAdd){
  getConverterMap().insert(ConverterPair(
    validator->getXMLTypeName(), converterToAdd));
}


RCP<const ValidatorXMLConverter>
ValidatorXMLConverterDB::getConverter(const ParameterEntryValidator& validator)
{
  ConverterMap::const_iterator it = getConverterMap().find(validator.getXMLTypeName());
  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindValidatorConverterException,
    "Could not find a ValidatorXMLConverter for validator type " <<
     validator.getXMLTypeName() << std::endl <<
     "Try adding an appropriate converter to the ValidatorXMLConverterDB " <<
     "in order solve this problem." << std::endl << std::endl
  )
  return it->second;
}


RCP<const ValidatorXMLConverter>
ValidatorXMLConverterDB::getConverter(const XMLObject& xmlObject)
{
  std::string validatorType = xmlObject.getRequired(
    ValidatorXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(validatorType);
  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindValidatorConverterException,
    "Could not find a ValidatorXMLConverter for type " << validatorType <<
    std::endl <<
    "Try adding an appropriate converter to the ValidatorXMLConverterDB " <<
    "in order solve this problem." << std::endl << std::endl
  )
  return it->second;
}


XMLObject ValidatorXMLConverterDB::convertValidator(
  RCP<const ParameterEntryValidator> validator,
  const ValidatortoIDMap& validatorIDsMap,
  bool assignID)
{
  return getConverter(*validator)->fromValidatortoXML(
    validator, validatorIDsMap, assignID);
}


RCP<ParameterEntryValidator> ValidatorXMLConverterDB::convertXML(
  const XMLObject& xmlObject,
  const IDtoValidatorMap& validatorIDsMap)
{
  return ValidatorXMLConverterDB::
    getConverter(xmlObject)->fromXMLtoValidator(xmlObject, validatorIDsMap);
}


ValidatorXMLConverterDB::ConverterMap&
ValidatorXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  return masterMap;
  // See default setup code below!
}


void ValidatorXMLConverterDB::printKnownConverters(std::ostream& out){
  out << "Known ValidatorXMLConverters: " << std::endl;
  for(
    ConverterMap::const_iterator it = getConverterMap().begin();
    it != getConverterMap().end();
    ++it)
  {
    out << "\t" << it->first <<std::endl;
  }
}


} // namespace Teuchos


namespace {


TEUCHOS_STATIC_SETUP()
{
  TEUCHOS_ADD_NUMBERTYPE_VALIDATOR_CONVERTERS(int);
  TEUCHOS_ADD_ENHANCEDNUMBERVALIDATOR_CONVERTER(double);
  TEUCHOS_ADD_ENHANCEDNUMBERVALIDATOR_CONVERTER(float);

  TEUCHOS_ADD_ARRAYVALIDATOR_CONVERTER(Teuchos::EnhancedNumberValidator<double>, double);
  TEUCHOS_ADD_ARRAYVALIDATOR_CONVERTER(Teuchos::EnhancedNumberValidator<float>, float);

  TEUCHOS_ADD_ARRAYVALIDATOR_CONVERTER(Teuchos::FileNameValidator, std::string);
  TEUCHOS_ADD_ARRAYVALIDATOR_CONVERTER(Teuchos::StringValidator, std::string);

  TEUCHOS_ADD_NUMBERTYPE_VALIDATOR_CONVERTERS(long long int);

  TEUCHOS_ADD_STRINGTOINTEGRALVALIDATOR_CONVERTER(Teuchos::EVerbosityLevel);

  TEUCHOS_ADD_VALIDATOR_CONVERTER(Teuchos::FileNameValidator, Teuchos::FileNameValidatorXMLConverter);
  TEUCHOS_ADD_VALIDATOR_CONVERTER(Teuchos::StringValidator, Teuchos::StringValidatorXMLConverter);
  TEUCHOS_ADD_VALIDATOR_CONVERTER(Teuchos::AnyNumberParameterEntryValidator, Teuchos::AnyNumberValidatorXMLConverter);
  TEUCHOS_ADD_VALIDATOR_CONVERTER(Teuchos::BoolParameterEntryValidator, Teuchos::BoolValidatorXMLConverter);

}


} // namespace
