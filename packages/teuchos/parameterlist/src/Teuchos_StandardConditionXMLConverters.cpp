// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_StandardConditionXMLConverters.hpp"
#include "Teuchos_ConditionXMLConverterDB.hpp"
#include "Teuchos_XMLConditionExceptions.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"

namespace Teuchos{

RCP<Condition> BoolLogicConditionConverter::convertXML(
    const XMLObject& xmlObj,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const
{
  Condition::ConstConditionList conditions;
  for(int i = 0; i < xmlObj.numChildren(); ++i){
    conditions.push_back(
      ConditionXMLConverterDB::convertXML(xmlObj.getChild(i), entryIDsMap));
  }
  return getSpecificBoolLogicCondition(conditions);
}

void BoolLogicConditionConverter::convertCondition(
  const RCP<const Condition> condition,
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const
{
  RCP<const BoolLogicCondition> castedCon =
    rcp_dynamic_cast<const BoolLogicCondition>(condition, true);

  const Condition::ConstConditionList conditions = castedCon->getConditions();
  for(
    Condition::ConstConditionList::const_iterator it = conditions.begin();
    it != conditions.end();
    ++it)
  {
    xmlObj.addChild(ConditionXMLConverterDB::convertCondition(*it, entryIDsMap));
  }
}

RCP<BoolLogicCondition>
OrConditionConverter::getSpecificBoolLogicCondition(
  Condition::ConstConditionList& conditions) const
{
  return rcp( new OrCondition(conditions));
}

RCP<BoolLogicCondition>
AndConditionConverter::getSpecificBoolLogicCondition(
  Condition::ConstConditionList& conditions) const
{
  return rcp( new AndCondition(conditions));
}

RCP<BoolLogicCondition>
EqualsConditionConverter::getSpecificBoolLogicCondition(
  Condition::ConstConditionList& conditions) const
{
  return rcp( new EqualsCondition(conditions));
}

RCP<Condition> NotConditionConverter::convertXML(
    const XMLObject& xmlObj,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const
{
  return rcp(new NotCondition(
    ConditionXMLConverterDB::convertXML(xmlObj.getChild(0), entryIDsMap)));
}

void NotConditionConverter::convertCondition(
  const RCP<const Condition> condition,
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const
{
  RCP<const NotCondition> castedCondition =
    rcp_dynamic_cast<const NotCondition>(condition);
  xmlObj.addChild(ConditionXMLConverterDB::convertCondition(
    castedCondition->getChildCondition(), entryIDsMap));
}

RCP<Condition> ParameterConditionConverter::convertXML(
    const XMLObject& xmlObj,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const
{
  ParameterEntry::ParameterEntryID paramID =
    xmlObj.getRequired<ParameterEntry::ParameterEntryID>(
      getParameterEntryIdAttributeName());
  TEUCHOS_TEST_FOR_EXCEPTION(
    entryIDsMap.find(paramID) == entryIDsMap.end(),
    MissingParameterEntryDefinitionException,
    "Can't find a parameter entry with id " << paramID << " in the "
    "given entryIDsMap!" << std::endl << std::endl);
  return getSpecificParameterCondition(
    xmlObj, entryIDsMap.find(paramID)->second);
}

void ParameterConditionConverter::convertCondition(
  const RCP<const Condition> condition,
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const
{
  RCP<const ParameterCondition> castedCondition =
    rcp_dynamic_cast<const ParameterCondition>(condition, true);

  TEUCHOS_TEST_FOR_EXCEPTION(
    entryIDsMap.find(castedCondition->getParameter()) == entryIDsMap.end(),
    MissingParameterEntryDefinitionException,
    "Couldn't find an id for the parameter in the given entryIDsMap!" <<
    std::endl << std::endl);

  xmlObj.addAttribute(
    getParameterEntryIdAttributeName(),
    entryIDsMap.find(castedCondition->getParameter())->second);

  addSpecificXMLTraits(castedCondition, xmlObj);
}

RCP<ParameterCondition>
StringConditionConverter::getSpecificParameterCondition(
  const XMLObject& xmlObj,
  RCP<ParameterEntry> parameterEntry) const
{
  StringCondition::ValueList values;
  int result = xmlObj.findFirstChild(getValuesTagName());
  TEUCHOS_TEST_FOR_EXCEPTION(result == -1,
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
  return rcp(new StringCondition(parameterEntry, values));
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
  xmlObj.addChild(valueTag);
}

RCP<ParameterCondition>
BoolConditionConverter::getSpecificParameterCondition(
  const XMLObject& /* xmlObj */,
  RCP<ParameterEntry> parameterEntry) const
{
  return rcp(new BoolCondition(parameterEntry));
}

void BoolConditionConverter::addSpecificXMLTraits(
  RCP<const ParameterCondition> /* condition */, XMLObject& /* xmlObj */) const
{}


} //namespace Teuchos


