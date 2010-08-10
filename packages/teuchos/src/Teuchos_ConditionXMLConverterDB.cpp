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
#include "Teuchos_StandardParameterEntryConditions.hpp"



namespace Teuchos {


#define ADD_NUMBERCONDITION(T) \
  \
  NumberCondition< T > ##T##NumberCondition; \
  masterMap.insert( \
    ConverterPair(##T##NumberCondition.getTypeAttributeValue(), \
      rcp(new NumberConditionXMLConverter< T >)));


void ConditionXMLConverterDB::addConverter(
  ParameterEntryValidator& condition,
  RCP<ConditionXMLConverter> converterToAdd){
  getConverterMap().insert(
    ConverterPair(condition.getTypeAttributeValue(), converterToAdd));
}


RCP<const ConditionXMLConverter>
ConditionXMLConverterDB::getConverter(const Condition& condition){
  ConverterMap::const_iterator it = 
    getConverterMap().find(condition.getTypeAttributeValue());
  TEST_FOR_EXCEPTION(it != getConverterMap().end(),
    CantFindConditionConverterException,
    "Could not find a ConditionXMLConverter for a condition"
  )
  return it->second;
}


RCP<const ConditionXMLConverter>
ConditionXMLConverterDB::getConverter(const XMLObject& xmlObject)
{ 
  std::string conditionType = xmlObject.getRequired(
    ConditionXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(conditionType);
  TEST_FOR_EXCEPTION(it != getConverterMap().end(),
    CantFindConditionConverterException,
    "Could not find a ConditionXMLConverter for a condition"
  )
  return it->second;
}

XMLObject ConditionXMLConverterDB::convertCondition(
  RCP<const Condition> condition)
{
  return getConverter(*condition)->fromConditiontoXML(condition);
}
 
RCP<ParameterEntryCondition> ConditionXMLConverterDB::convertXML(
  const XMLObject& xmlObject)
{
  return ConditionXMLConverterDB::getConverter(xmlObject)->
    fromXMLtoCondition(xmlObject);
}

ConditionXMLConverterDB::ConverterMap&
ConditionXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  if(masterMap.size() == 0){

    ADD_NUMBERCONVERTER(int);
    typedef unsigned int uint;
    ADD_NUMBERCONVERTER(uint);
    typedef short int sint;
    ADD_NUMBERCONVERTER(sint);
    typedef unsigned short int usint;
    ADD_NUMBERCONVERTER(usint);
    typedef long int lint;
    ADD_NUMBERCONVERTER(lint);
    typedef unsigned long int ulint;
    ADD_NUMBERCONVERTER(ulint);
    ADD_NUMBERCONVERTER(double);
    ADD_NUMBERCONVERTER(float);

    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    typedef long long int llint;
    ADD_NUMBERCONVERTER(llint);
    typedef unsigned long long int ullint;
    ADD_NUMBERCONVERTER(ullint);
    #endif // HAVE_TEUCHOS_LONG_LONG_INT

    StringCondition stringCondition;
    masterMap.insert(
      ConverterPair(stringCondition.getTypeAttributeValue(), 
        rcp(new StringConditionXMLConverter)));

    BoolCondition boolCondition;
    masterMap.insert(
      ConverterPair(boolCondition.getTypeAttributeValue(), 
        rcp(new BoolConditionXMLConverter)));

    OrCondition orCondition;
    masterMap.insert(
      ConverterPair(orCondition.getTypeAttributeValue(), 
        rcp(new OrConditionXMLConverter)));

    AndCondition andCondition;
    masterMap.insert(
      ConverterPair(andCondition.getTypeAttributeValue(), 
        rcp(new AndConditionXMLConverter)));

    EqualsCondition equalsCondition;
    masterMap.insert(
      ConverterPair(equalsCondition.getTypeAttributeValue(), 
        rcp(new EqualsConditionXMLConverter)));

    NotCondition notCondition;
    masterMap.insert(
      ConverterPair(notCondition.getTypeAttributeValue(), 
        rcp(new NotConditionXMLConverter)));

  }
  return masterMap;
}


} // namespace Teuchos
