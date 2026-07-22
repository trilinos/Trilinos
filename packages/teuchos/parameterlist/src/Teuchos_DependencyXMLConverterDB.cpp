// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_DependencyXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"



namespace Teuchos {



void DependencyXMLConverterDB::addConverter(
  RCP<const Dependency> dependency,
  RCP<DependencyXMLConverter> converterToAdd)
{
  getConverterMap().insert(
    ConverterPair(dependency->getTypeAttributeValue(), converterToAdd));
}


RCP<const DependencyXMLConverter>
DependencyXMLConverterDB::getConverter(const Dependency& dependency)
{
  ConverterMap::const_iterator it =
    getConverterMap().find(dependency.getTypeAttributeValue());
  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindDependencyConverterException,
    "Could not find a DependencyXMLConverter for a dependency with "
    "attribute tag " << dependency.getTypeAttributeValue() <<
    "!" << std::endl <<
    "Try adding an appropriate converter to the DependencyXMLConverterDB " <<
    "in order to solve this problem." << std::endl << std::endl);
  return it->second;
}


RCP<const DependencyXMLConverter>
DependencyXMLConverterDB::getConverter(const XMLObject& xmlObject)
{
  std::string dependencyType = xmlObject.getRequired(
    DependencyXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(dependencyType);
  #ifdef HAVE_TEUCHOS_DEBUG
  std::ostringstream sout;
  printKnownConverters(sout);
  #endif
  std::string extraError =
  #ifdef HAVE_TEUCHOS_DEBUG
  sout.str();
  #else
  "";
  #endif

  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindDependencyConverterException,
    "Could not find a DependencyXMLConverter for a dependency of type " <<
    dependencyType << "!" << std::endl <<
    "Try adding an appropriate converter to the DependencyXMLConverterDB " <<
    "in order to solve this problem." << std::endl << std::endl << extraError
    );
  return it->second;
}

XMLObject DependencyXMLConverterDB::convertDependency(
  RCP<const Dependency> dependency,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
  ValidatortoIDMap& validatorIDsMap)
{
  return getConverter(*dependency)->fromDependencytoXML(
    dependency, entryIDsMap, validatorIDsMap);
}

RCP<Dependency> DependencyXMLConverterDB::convertXML(
    const XMLObject& xmlObject,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap)
{
  return DependencyXMLConverterDB::getConverter(xmlObject)->
    fromXMLtoDependency(xmlObject, entryIDsMap, validatorIDsMap);
}

DependencyXMLConverterDB::ConverterMap&
DependencyXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  return masterMap;
}


} // namespace Teuchos


namespace {


TEUCHOS_STATIC_SETUP()
{
    TEUCHOS_ADD_TEMPLATED_NUMBER_DEPS(int)
    TEUCHOS_ADD_NUMBER_VISUAL_DEP(float)
    TEUCHOS_ADD_RANGE_VALIDATOR_DEP(float)
    TEUCHOS_ADD_NUMBER_VISUAL_DEP(double)
    TEUCHOS_ADD_RANGE_VALIDATOR_DEP(double)

    TEUCHOS_ADD_TEMPLATED_NUMBER_DEPS(long long int)

    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::StringValidatorDependency,
      Teuchos::StringValidatorDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::StringVisualDependency,
      Teuchos::StringVisualDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::BoolValidatorDependency,
      Teuchos::BoolValidatorDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::BoolVisualDependency,
      Teuchos::BoolVisualDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::ConditionVisualDependency,
      Teuchos::ConditionVisualDependencyXMLConverter)
}


} //namespace

