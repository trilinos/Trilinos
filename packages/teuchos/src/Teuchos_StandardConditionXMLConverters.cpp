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
#include "Teuchos_ConditionXMLConveterDB.hpp"

namespace Teuchos{

RCP<Condition> BinaryLogicalConditionConverter::convertXML(
    const XMLObject& xmlObj) const
{
  Condtion::CondtionList conditions;
  for(int i = 0; i < xmlObj.numChildren(); ++i){
    conditions.push_back(ConditionXMLConverterDB::convertXML(xmlObj.getChild(i)));
  }
  return getSpecificBinaryLogicalCondition(conditions); 
}

void BinaryLogicalConditionConverter::convertCondition(
  const RCP<const Condition> condition, XMLObject& xmlObj) const
{
  RCP<const BinaryLogicalCondition> castedCon = 
    rcp_dynamic_cast<const BinaryLogicalCondition>(condition, true);

  const Condition::ConditionList conditions = castedCon->getConditions();
  for(
    Condition::ConditionList::const_iterator it = conditions.begin();
    it != conditions.end();
    ++it)
  {
    xmlObj.addChild(ConditionXMLConverterDB::convertCondition(*it));    
  }
}

RCP<BinaryLogicalCondition> 
OrConditionConverter::getSpecificBinaryLogicalCondition(
  Condition::ConditionList& conditions) const
{
  return rcp( new OrCondition(conditions));    
}

RCP<BinaryLogicalCondition> 
AndConditionConverter::getSpecificBinaryLogicalCondition(
  Condition::ConditionList& conditions) const
{
  return rcp( new AndCondition(conditions));    
}

RCP<BinaryLogicalCondition> 
EqualsConditionConverter::getSpecificBinaryLogicalCondition(
  Condition::ConditionList& conditions) const
{
  return rcp( new EqualsCondition(conditions));    
}

RCP<Condition> NotConditionConverter::convertXML(const XMLObject& xmlObj) const{
  return rcp(new NotCondition(
    ConditionXMLConverterDB::convertXML(xmlObj.getChild(0))));
}

void NotConditionConverter::convertCondition(
  const RCP<const Condition> condition,
  XMLObject& xmlObj) const
{
  RCP<const NotCondition> castedCondition =
    rcp_dynamic_cast<const NotCondition>(condition);
  XMLObject.addChild(ConditionXMLConverterDB::convertCondition(
    castedCondition->getChildCondition()));
}

RCP<Condition> ParameterConditionConverter::convertXML(
  const XMLObject& xmlObj) const
{
   std::string parameterName = xmlObj.getRequired(getParameterNameAttributeName());
   RCP<ParameterList> parentList
}

void ParameterConditionConverter::convertCondition(
  const RCP<const Condition> condition, 
  XMLObject& xmlObj) const
{
  RCP<const ParameterCondtion> castedCondition = 
    rcp_dynamic_cast<const ParameterCondtion>(condition, true);
  xmlObj.addAttribute(getParameterNameAttributeName(), castedCondition->getParameterName());
  xmlObj.addAttribute(getParentListNameAttributeName(), castedCondition->getParentList()->name());
  xmlObj.addBool(getWhenParamEqualsValueAttributeName(), castedCondition->getWhenParamEqualsValue());
  
  addSpecificXMLTraits(castedCondition, xmlObj);
}

} //namespace Teuchos


