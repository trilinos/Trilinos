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
  const XMLObject& xmlObj,
  RCP<ParameterEntry> parameterEntry) const
{
  return rcp(new BoolCondition(parameterEntry));
}

void BoolConditionConverter::addSpecificXMLTraits(
  RCP<const ParameterCondition> condition, XMLObject& xmlObj) const             
{}
 

} //namespace Teuchos


