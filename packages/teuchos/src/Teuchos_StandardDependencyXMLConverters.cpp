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

#include "Teuchos_StandardDependencyXMLConverters.hpp"


namespace Teuchos{


RCP<Dependency> VisualDependencyConverter::convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents) const
{
  bool showIf = xmlObj.getWithDefault(getShowIfAttributeName, true);
  return convertSpecialVisualAttributes(xmlObj,
    dependees,
    dependents,
    showIf);
}

void VisualDependencyConverter::convertDependency(
  const RCP<const Dependency> dependency, 
  XMLObject& xmlObj) const
{
  RCP<const VisualDependency> castedDep = 
    rcp_dynamic_cast<const VisualDependency>(dependency, true);

  xmlObj.addBool(getShowIfAttributeName(), castedDep->getShowIf());
  convertSpecialVisualAttributes(dependency, xmlObj);
}

RCP<Dependency> ValidatorDependencyXMLConverter::convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A Validator Dependency can only have 1 dependee!" << 
    std::endl << std::endl);
  return convertSpecialValidatorAttributes(
    xmlObj, *(dependees.begin()), dependents);
}

void ValidatorDependencyXMLConverter::convertDependency(
  const RCP<const Dependency> dependency, 
  XMLObject& xmlObj) const
{
  convertSpecialValidatorAttributes(dependency, xmlObj);
}
  
 
StringVisualDependencyConverter::convertSpecialVisualAttributes(
  RCP<const VisualDepenency> dependency, XMLObject& xmlObj) const
{
  RCP<const StringVisualDependency> castedDependency = 
    rcp_dynamic_cast<const StringVisualDependency(dependency, true);
  StringVisualDependency::ValueList valueList = castedDependency->getValues();
  XMLObject valuesTag(getStringValuesTagName());
  for(ValueList::const_interator it = valueList.begin(); it != valueList.end(); ++it){
    XMLObject stringValue(getStringTagName());
    stringValue.addAttribute(getValueAttributeName(), *it);
    valuesTag.addChild(stringValue);
  }
  xmlObj.addChild(valuesTag);
}
  
RCP<VisualDepenency> StringVisualDependencyConverter::convertSpecialVisualAttributes(
  XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A StringVisualDependency can only have 1 dependee!" << std::endl << std::endl
  );

  StringVisualDependency::ValueList valueList;
  XMLObject valuesTag = xmlObj.findFirstChild(getStringValuesTagName());
  
  TEST_FOR_EXCEPTION(nonnull(valuesTag),
    ValuesTagMissingException,
    "Couldn't find " << getStringValuesTagName() << " tag for a " <<
    "StringVisualDependency!" << std::endl <<std::endl);

  for(int i=0; i<valuesTag.numChildren(); ++i){
    XMLObject child = valuesTag.getChild(i);
    valueList.push_back(child.getRequired(getValueAttributeName()));
  }

  return rcp(new StringVisualDependency(dependees, dependents, valueList, showIf));
}
  
void BoolVisualDependencyConverter::convertSpecialVisualAttributes(
  RCP<const VisualDepenency> dependency,
  XMLObject& xmlObj) const {}
  
RCP<VisualDepenency> 
BoolVisualDependencyConverter::convertSpecialVisualAttributes(
  XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A StringVisualDependency can only have 1 dependee!" <<
    std::endl << std::endl);
  return rcp(new BoolVisualDependency(
    *(dependees.begin()), dependents, showIf));
}

void ConditionVisualDependencyConverter::convertSpecialVisualAttributes(
  RCP<const VisualDepenency> dependency,
  XMLObject& xmlObj) const
{
  RCP<ConditionVisualDependency> castedDependency = 
    rcp_dynamic_cast<ConditionVisualDependency>(dependency);
  xmlObj.addChild(
    ConditionXMLConvertDB::convertCondition(castedDependency->getCondition()));
}
  
RCP<VisualDepenency> 
ConditionVisualDependency::convertSpecialVisualAttributes(
  XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf) const
{
  XMLObject conditionObj = xmlObj.findFirstChild(Condition::getXMLTagName());
  TEST_FOR_EXCEPTION(conditionObj == NULL,
    MissingConditionTagException,
    "ConditionVisualDependencies must have a Condition tag!"
  );
  Teuchos::RCP<Condition> condition = 
    ConditionXMLConvertDB::convertXML(conditionObj);
  return rcp(new ConditionVisualDependency(
    ConditionXMLConvertDB::convertXML(condition, dependents, showIf)));
}
 
void
StringValidatorDependencyConverter::convertSpecialValidatorAttributes(
  RCP<const ValidatorDepenency> dependency,
  XMLObject& xmlObj) const
{
  RCP<const StringValidatorDependency> castedDependency = 
    rcp_dynamic_cast<const StringValidatorDependency>(dependency, true);  
  XMLObject valueMapTag(getValuesAndValidatorsTag()); 
  const StringVisualDependency::ValueToValidatorMap valuesAndValidators = 
    castedDependency->getValuesAndValidators();
  for(
    StringValidatorDependencyConverter::ValidatorDepenency::const_iterator it =
      valuesAndValidators.begin();
    it != valuesAndValidators.end();
    ++it)
  {
    XMLObject pairTag(getPairTag());
    pairTag.addAttribute(getValueAttributeName(), it->first);
    pairTag.addAttribute(getValidatorIDAttributeName(), 
      ParameterEntryValidator::getValidatorID(it->second));
    valueMapTag.addChild(pairTag);
  }  
}

RCP<ValidatorDepenency> 
StringValidatorDependencyConverter::convertSpecialValidatorAttributes(
  XMLObject& xmlObj,
  const RCP<const ParameterEntry> dependee,
  const Dependency::ParameterEntryList dependents) const
{
  ValueToValidatorMap valueValidatorMap;
  XMLObject valuesAndValidatorTag = 
    xmlObj.findFirstChild(getValuesAndValidatorsTag());
  TEST_FOR_EXCEPTION(valuesAndValidatorTag != NULL,
    MissingValuesAndValidatorsTagException,
    "Error: All StringValidatorDependencies must have a " << 
    getValuesAndValidatorsTag() << "tag!" << std::endl << std::endl);
  for(int i=0; i < valuesAndValidatorTag.numChildren(); ++i){
    XMLObject child = valuesAndValidatorTag.getChild(i);
    std::string value = child.getRequired(getValueAttributeName());
    RCP<ParameterEntryValidator> validator = 
      ParameterEntryValidator::getValidator(child.getRequired(
        getValidatorIDAttributeName()));
    valueValidatorMap.insert(
      StringValidatorDependency::ValueToValidatorPair(value, validator));
  }

  RCP<ParameterEntryValidator> defaultValidator = null;
  if(xmlObj.hasAttribute(getDefaultValidatorIDAttributeName())){
    defaultValidator = ParameterEntryValidator::getValidator(
      xmlObj.getRequired(getDefaultValidatorIDAttributeName()));
  }

  return rcp(new StringValidatorDependency(
    dependee, dependents, valueValidatorMap, defaultValidator));
}

void
BoolValidatorDependencyConverter::convertSpecialValidatorAttributes(
  RCP<const ValidatorDepenency> dependency,
  XMLObject& xmlObj) const
{
  RCP<const BoolValidatorDependency> castedDependency = 
    rcp_dynamic_cast<const BoolValidatorDependency>(dependency, true);  
  xmlObj.addAttribute(getFalseValidatorIDAttributeName(),
    ParameterEntryValidator::getValidatorID(
      castedDependency->getFalseValidator()));

  xmlObj.addAttribute(getTrueValidatorIDAttributeName(),
    ParameterEntryValidator::getValidatorID(
      castedDependency->getTrueValidator()));
}

RCP<ValidatorDepenency> 
BoolValidatorDependencyConverter::convertSpecialValidatorAttributes(
  XMLObject& xmlObj,
  const RCP<const ParameterEntry> dependee,
  const Dependency::ParameterEntryList dependents) const
{
  RCP<ParameterEntryValidator> falseValidator = 
    ParameterEntryValidator::getValidator(
      xmlObj.getRequired(getFalseValidatorIDAttributeName()));
    
  RCP<ParameterEntryValidator> trueValidator = 
    ParameterEntryValidator::getValidator(
      xmlObj.getRequired(getTrueValidatorIDAttributeName()));

  return rcp(new BoolValidatorDependency(
    dependee, dependents, falseValidator, trueValidator));
}



} //namespace Teuchos

