// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
