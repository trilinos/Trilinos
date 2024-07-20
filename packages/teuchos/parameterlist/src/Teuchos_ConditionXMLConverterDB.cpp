// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConditionXMLConverterDB.hpp"
#include "Teuchos_XMLConditionExceptions.hpp"
#include "Teuchos_StaticSetupMacro.hpp"



namespace Teuchos {


void ConditionXMLConverterDB::addConverter(
  RCP<const Condition> condition,
  RCP<ConditionXMLConverter> converterToAdd){
  getConverterMap().insert(
    ConverterPair(condition->getTypeAttributeValue(), converterToAdd));
}


RCP<const ConditionXMLConverter>
ConditionXMLConverterDB::getConverter(const Condition& condition){
  ConverterMap::const_iterator it =
    getConverterMap().find(condition.getTypeAttributeValue());
  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindConditionConverterException,
    "Could not find a ConditionXMLConverter for a condition of type " <<
    condition.getTypeAttributeValue() << " when writing out a condition to " <<
    "xml." << std::endl << std::endl
  )
  return it->second;
}


RCP<const ConditionXMLConverter>
ConditionXMLConverterDB::getConverter(const XMLObject& xmlObject)
{
  std::string conditionType = xmlObject.getRequired(
    ConditionXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(conditionType);
  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindConditionConverterException,
    "Could not find a ConditionXMLConverter for a condition of type " <<
    conditionType << " when reading in a condition from " <<
    "xml." << std::endl << std::endl
  )
  return it->second;
}

XMLObject ConditionXMLConverterDB::convertCondition(
  RCP<const Condition> condition,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap)
{
  return getConverter(*condition)->fromConditiontoXML(condition, entryIDsMap);
}

RCP<Condition> ConditionXMLConverterDB::convertXML(
  const XMLObject& xmlObject,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap)
{
  return ConditionXMLConverterDB::getConverter(xmlObject)->
    fromXMLtoCondition(xmlObject, entryIDsMap);
}

ConditionXMLConverterDB::ConverterMap&
ConditionXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  return masterMap;
}


} // namespace Teuchos


namespace {


TEUCHOS_STATIC_SETUP()
{
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(int);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(unsigned int);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(short int);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(unsigned short int);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(long int);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(unsigned long int);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(double);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(float);

    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(long long int);
    TEUCHOS_ADD_NUMBERCONDITION_CONVERTER(unsigned long long int);

    Teuchos::ConditionXMLConverterDB::addConverter(
      Teuchos::DummyObjectGetter<Teuchos::StringCondition>::
        getDummyObject(),
      Teuchos::rcp(new Teuchos::StringConditionConverter));

    Teuchos::ConditionXMLConverterDB::addConverter(
      Teuchos::DummyObjectGetter<Teuchos::BoolCondition>::
        getDummyObject(),
      Teuchos::rcp(new Teuchos::BoolConditionConverter));

    Teuchos::ConditionXMLConverterDB::addConverter(
      Teuchos::DummyObjectGetter<Teuchos::OrCondition>::
        getDummyObject(),
      Teuchos::rcp(new Teuchos::OrConditionConverter));

    Teuchos::ConditionXMLConverterDB::addConverter(
      Teuchos::DummyObjectGetter<Teuchos::AndCondition>::
        getDummyObject(),
      Teuchos::rcp(new Teuchos::AndConditionConverter));

    Teuchos::ConditionXMLConverterDB::addConverter(
      Teuchos::DummyObjectGetter<Teuchos::EqualsCondition>::
        getDummyObject(),
      Teuchos::rcp(new Teuchos::EqualsConditionConverter));

    Teuchos::ConditionXMLConverterDB::addConverter(
      Teuchos::DummyObjectGetter<Teuchos::NotCondition>::
        getDummyObject(),
      Teuchos::rcp(new Teuchos::NotConditionConverter));
}


} // namespace
