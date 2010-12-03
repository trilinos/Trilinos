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

#include "Teuchos_ConditionXMLConverterDB.hpp"
#include "Teuchos_StandardConditionXMLConverters.hpp"
#include "Teuchos_StandardConditions.hpp"
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
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
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
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
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
    TEUCHOS_ADD_NUMBERCONVERTER(int);
    TEUCHOS_ADD_NUMBERCONVERTER(unsigned int);
    TEUCHOS_ADD_NUMBERCONVERTER(short int);
    TEUCHOS_ADD_NUMBERCONVERTER(unsigned short int);
    TEUCHOS_ADD_NUMBERCONVERTER(long int);
    TEUCHOS_ADD_NUMBERCONVERTER(unsigned long int);
    TEUCHOS_ADD_NUMBERCONVERTER(double);
    TEUCHOS_ADD_NUMBERCONVERTER(float);

    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    TEUCHOS_ADD_NUMBERCONVERTER(long long int);
    TEUCHOS_ADD_NUMBERCONVERTER(unsigned long long int);
    #endif // HAVE_TEUCHOS_LONG_LONG_INT

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
