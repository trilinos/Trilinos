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

#include "Teuchos_StandardConditionXMLConverters.hpp"
#include "Teuchos_ConditionXMLConverterDB.hpp"
#include "Teuchos_XMLConditionExceptions.hpp"

namespace Teuchos{

RCP<Condition> BinaryLogicalConditionConverter::convertXML(
    const XMLObject& xmlObj) const
{
  Condition::ConstConditionList conditions;
  for(int i = 0; i < xmlObj.numChildren(); ++i){
    conditions.push_back(
      ConditionXMLConverterDB::convertXML(xmlObj.getChild(i)));
  }
  return getSpecificBinaryLogicalCondition(conditions); 
}

void BinaryLogicalConditionConverter::convertCondition(
  const RCP<const Condition> condition, XMLObject& xmlObj) const
{
  RCP<const BinaryLogicalCondition> castedCon = 
    rcp_dynamic_cast<const BinaryLogicalCondition>(condition, true);

  const Condition::ConstConditionList conditions = castedCon->getConditions();
  for(
    Condition::ConstConditionList::const_iterator it = conditions.begin();
    it != conditions.end();
    ++it)
  {
    xmlObj.addChild(ConditionXMLConverterDB::convertCondition(*it));    
  }
}

RCP<BinaryLogicalCondition> 
OrConditionConverter::getSpecificBinaryLogicalCondition(
  Condition::ConstConditionList& conditions) const
{
  return rcp( new OrCondition(conditions));    
}

RCP<BinaryLogicalCondition> 
AndConditionConverter::getSpecificBinaryLogicalCondition(
  Condition::ConstConditionList& conditions) const
{
  return rcp( new AndCondition(conditions));    
}

RCP<BinaryLogicalCondition> 
EqualsConditionConverter::getSpecificBinaryLogicalCondition(
  Condition::ConstConditionList& conditions) const
{
  return rcp( new EqualsCondition(conditions));    
}

RCP<Condition> NotConditionConverter::convertXML(
  const XMLObject& xmlObj) const
{
  return rcp(new NotCondition(
    ConditionXMLConverterDB::convertXML(xmlObj.getChild(0))));
}

void NotConditionConverter::convertCondition(
  const RCP<const Condition> condition,
  XMLObject& xmlObj) const
{
  RCP<const NotCondition> castedCondition =
    rcp_dynamic_cast<const NotCondition>(condition);
  xmlObj.addChild(ConditionXMLConverterDB::convertCondition(
    castedCondition->getChildCondition()));
}

RCP<Condition> ParameterConditionConverter::convertXML(
  const XMLObject& xmlObj) const
{
  ParameterEntry::ParameterEntryID paramID = 
    xmlObj.getRequired<ParameterEntry::ParameterEntryID>(
      getParameterEntryIDAttributeName());
  bool whenParamEqualsValue = 
    xmlObj.getRequiredBool(getWhenParamEqualsValueAttributeName());
  return getSpecificParameterCondition(
    xmlObj, paramID, whenParamEqualsValue);
}

void ParameterConditionConverter::convertCondition(
  const RCP<const Condition> condition, 
  XMLObject& xmlObj) const
{
  RCP<const ParameterCondition> castedCondition = 
    rcp_dynamic_cast<const ParameterCondition>(condition, true);

  xmlObj.addAttribute(
    getParameterEntryIDAttributeName(),
    ParameterEntry::getParameterEntryID(castedCondition->getParameter()));

  xmlObj.addBool(
    getWhenParamEqualsValueAttributeName(), 
    castedCondition->getWhenParamEqualsValue());
  
  addSpecificXMLTraits(castedCondition, xmlObj);
}

RCP<ParameterCondition> 
StringConditionConverter::getSpecificParameterCondition(
  const XMLObject& xmlObj,
  ParameterEntry::ParameterEntryID paramID,
  bool whenParamEqualsValue) const
{
  StringCondition::ValueList values;
  int result = xmlObj.findFirstChild(getValuesTagName());
  TEST_FOR_EXCEPTION(result == -1,
    MissingValuesTagException,
    "A StringCondtion must have a tag with the name " <<
    getValuesTagName() << " as one of it's children!");

  XMLObject valuesTag = xmlObj.getChild(result);
  for(int i=0; i< valuesTag.numChildren(); ++i){
    XMLObject child = valuesTag.getChild(i);
    if(child.getTag() == getStringTagName()){
      values.append(child.getRequired(getStringValueAttributeName()));
    }
  }
  return rcp(new StringCondition(
    ParameterEntry::getParameterEntry(paramID), values, whenParamEqualsValue));
}

void StringConditionConverter::addSpecificXMLTraits(
  RCP<const ParameterCondition> condition, XMLObject& xmlObj) const             
{
  RCP<const StringCondition> castedCon =
    rcp_dynamic_cast<const StringCondition>(condition, true);
  XMLObject valueTag(getValuesTagName());
  for(
    StringCondition::ValueList::const_iterator it =
      castedCon->getValueList().begin();
    it != castedCon->getValueList().end();
    ++it)
  {
    XMLObject stringTag(getStringTagName());
    stringTag.addAttribute(getStringValueAttributeName(), *it);
    valueTag.addChild(stringTag);
  }
}
 
RCP<ParameterCondition> 
BoolConditionConverter::getSpecificParameterCondition(
  const XMLObject& xmlObj,
  ParameterEntry::ParameterEntryID paramID,
  bool whenParamEqualsValue) const
{
  return rcp(new BoolCondition(
    ParameterEntry::getParameterEntry(paramID), whenParamEqualsValue));
}

void BoolConditionConverter::addSpecificXMLTraits(
  RCP<const ParameterCondition> condition, XMLObject& xmlObj) const             
{}
 

} //namespace Teuchos


