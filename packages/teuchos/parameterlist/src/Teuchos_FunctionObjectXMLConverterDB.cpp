// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_FunctionObjectXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"
#include "Teuchos_XMLFunctionObjectExceptions.hpp"



namespace Teuchos {


void FunctionObjectXMLConverterDB::addConverter(
  RCP<const FunctionObject> function,
  RCP<FunctionObjectXMLConverter> converterToAdd){
  getConverterMap().insert(
    ConverterPair(function->getTypeAttributeValue(), converterToAdd));
}


RCP<const FunctionObjectXMLConverter>
FunctionObjectXMLConverterDB::getConverter(const FunctionObject& function){
  ConverterMap::const_iterator it =
    getConverterMap().find(function.getTypeAttributeValue());
  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindFunctionObjectConverterException,
    "Could not find a FunctionObjectXMLConverter for a FuncitonObject of type " <<
    function.getTypeAttributeValue() << " when writing out a condition to " <<
    "xml." << std::endl << std::endl
  )
  return it->second;
}


RCP<const FunctionObjectXMLConverter>
FunctionObjectXMLConverterDB::getConverter(const XMLObject& xmlObject)
{
  std::string functionType = xmlObject.getRequired(
    FunctionObjectXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(functionType);
  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindFunctionObjectConverterException,
    "Could not find a FunctionObjectXMLConverter for a condition of type " <<
    functionType << " when reading in a condition from " <<
    "xml." << std::endl << std::endl
  )
  return it->second;
}

XMLObject FunctionObjectXMLConverterDB::convertFunctionObject(
  RCP<const FunctionObject> function)
{
  return getConverter(*function)->fromFunctionObjecttoXML(function);
}

RCP<FunctionObject> FunctionObjectXMLConverterDB::convertXML(
  const XMLObject& xmlObject)
{
  return FunctionObjectXMLConverterDB::getConverter(xmlObject)->
    fromXMLtoFunctionObject(xmlObject);
}

FunctionObjectXMLConverterDB::ConverterMap&
FunctionObjectXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  return masterMap;
}


} // namespace Teuchos


namespace {


TEUCHOS_STATIC_SETUP()
{

    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(short int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned short int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(long int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned long int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(double);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(float);

    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(long long int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned long long int);

}


} // namespace
