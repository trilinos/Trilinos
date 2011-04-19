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
#include "Teuchos_ConditionXMLConverterDB.hpp"


namespace Teuchos{


RCP<Dependency> VisualDependencyXMLConverter::convertXML(
  const XMLObject& xmlObj, 
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  bool showIf = xmlObj.getWithDefault(getShowIfAttributeName(), true);
  return convertSpecialVisualAttributes(
    xmlObj,
    dependees,
    dependents,
    showIf,
    entryIDsMap);
}

void VisualDependencyXMLConverter::convertDependency(
  const RCP<const Dependency> dependency, 
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
  ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const VisualDependency> castedDep = 
    rcp_dynamic_cast<const VisualDependency>(dependency, true);

  xmlObj.addBool(getShowIfAttributeName(), castedDep->getShowIf());
  convertSpecialVisualAttributes(castedDep, xmlObj, entryIDsMap);
}

RCP<Dependency> ValidatorDependencyXMLConverter::convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents,
    const XMLParameterListReader::EntryIDsMap& /*entryIDsMap*/,
    const IDtoValidatorMap& validatorIDsMap) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A Validator Dependency can only have 1 dependee!" << 
    std::endl << std::endl);
  return convertSpecialValidatorAttributes(
    xmlObj, *(dependees.begin()), dependents, validatorIDsMap);
}

void ValidatorDependencyXMLConverter::convertDependency(
  const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& /*entryIDsMap*/,
    ValidatortoIDMap& validatorIDsMap) const
{
  RCP<const ValidatorDependency> castedDep = 
    rcp_dynamic_cast<const ValidatorDependency>(dependency, true);
  convertSpecialValidatorAttributes(castedDep, xmlObj, validatorIDsMap);
}
  
 
void StringVisualDependencyXMLConverter::convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& /*entryIDsMap*/) const
{
  RCP<const StringVisualDependency> castedDependency = 
    rcp_dynamic_cast<const StringVisualDependency>(dependency, true);
  StringVisualDependency::ValueList valueList = castedDependency->getValues();
  XMLObject valuesTag(getStringValuesTagName());
  for(
    StringVisualDependency::ValueList::const_iterator it = valueList.begin(); 
    it != valueList.end(); 
    ++it)
  {
    XMLObject stringValue(getStringTagName());
    stringValue.addAttribute(getValueAttributeName(), *it);
    valuesTag.addChild(stringValue);
  }
  xmlObj.addChild(valuesTag);
}
  
RCP<VisualDependency> 
StringVisualDependencyXMLConverter::convertSpecialVisualAttributes(
  const XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf,
  const XMLParameterListReader::EntryIDsMap& /*entryIDsMap*/) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A StringVisualDependency can only have 1 dependee!" << 
    std::endl << std::endl);

  StringVisualDependency::ValueList valueList;
  int valuesTagIndex = xmlObj.findFirstChild(getStringValuesTagName());
  
  TEST_FOR_EXCEPTION(valuesTagIndex < 0,
    ValuesTagMissingException,
    "Couldn't find " << getStringValuesTagName() << " tag for a " <<
    "StringVisualDependency!" << std::endl <<std::endl);

  XMLObject valuesTag = xmlObj.getChild(valuesTagIndex);

  for(int i=0; i<valuesTag.numChildren(); ++i){
    XMLObject child = valuesTag.getChild(i);
    valueList.push_back(child.getRequired(getValueAttributeName()));
  }

  return rcp(
    new StringVisualDependency(
      *(dependees.begin()), 
      dependents, 
      valueList, 
      showIf));
}
  
void BoolVisualDependencyXMLConverter::convertSpecialVisualAttributes(
  RCP<const VisualDependency> dependency,
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& /*entryIDsMap*/) const
{}
  
RCP<VisualDependency> 
BoolVisualDependencyXMLConverter::convertSpecialVisualAttributes(
  const XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf,
  const XMLParameterListReader::EntryIDsMap& /*entryIDsMap*/) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A BoolVisualDependency can only have 1 dependee!" <<
    std::endl << std::endl);
  return rcp(new BoolVisualDependency(
    *(dependees.begin()), dependents, showIf));
}

void ConditionVisualDependencyXMLConverter::convertSpecialVisualAttributes(
  RCP<const VisualDependency> dependency,
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const
{
  RCP<const ConditionVisualDependency> castedDependency = 
    rcp_dynamic_cast<const ConditionVisualDependency>(dependency, true);
  xmlObj.addChild(
    ConditionXMLConverterDB::convertCondition(
      castedDependency->getCondition(), entryIDsMap));
}
  
RCP<VisualDependency> 
ConditionVisualDependencyXMLConverter::convertSpecialVisualAttributes(
  const XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap) const
{
  int conditionIndex = xmlObj.findFirstChild(Condition::getXMLTagName());
  TEST_FOR_EXCEPTION(conditionIndex < 0,
    MissingConditionTagException,
    "ConditionVisualDependencies must have a Condition tag!"
  );
  XMLObject conditionObj = xmlObj.getChild(conditionIndex);
  Teuchos::RCP<Condition> condition = 
    ConditionXMLConverterDB::convertXML(conditionObj, entryIDsMap);
  return rcp(new ConditionVisualDependency(condition, dependents, showIf));
}
 
void
StringValidatorDependencyXMLConverter::convertSpecialValidatorAttributes(
  RCP<const ValidatorDependency> dependency,
  XMLObject& xmlObj,
  ValidatortoIDMap& validatorIDsMap) const
{
  RCP<const StringValidatorDependency> castedDependency = 
    rcp_dynamic_cast<const StringValidatorDependency>(dependency, true);  
  XMLObject valueMapTag(getValuesAndValidatorsTag()); 
  const StringValidatorDependency::ValueToValidatorMap valuesAndValidators = 
    castedDependency->getValuesAndValidators();
  for(
    StringValidatorDependency::ValueToValidatorMap::const_iterator it =
      valuesAndValidators.begin();
    it != valuesAndValidators.end();
    ++it)
  {
    XMLObject pairTag(getPairTag());
    pairTag.addAttribute(getValueAttributeName(), it->first);
    if(validatorIDsMap.find(it->second) == validatorIDsMap.end()){
      validatorIDsMap.insert(it->second);
    }
    pairTag.addAttribute(getValidatorIdAttributeName(), 
      validatorIDsMap.find(it->second)->second);
    valueMapTag.addChild(pairTag);
  }  
  xmlObj.addChild(valueMapTag);
  RCP<const ParameterEntryValidator> defaultVali = 
    castedDependency->getDefaultValidator();
  if(nonnull(defaultVali)){
    if(validatorIDsMap.find(defaultVali) == validatorIDsMap.end()){
      validatorIDsMap.insert(defaultVali);
    }
    xmlObj.addAttribute(
      getDefaultValidatorIDAttributeName(),
      validatorIDsMap.find(defaultVali)->second);
  }
}

RCP<ValidatorDependency> 
StringValidatorDependencyXMLConverter::convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const
{
  StringValidatorDependency::ValueToValidatorMap valueValidatorMap;
  int valuesAndValidatorIndex = 
    xmlObj.findFirstChild(getValuesAndValidatorsTag()); 
    
  TEST_FOR_EXCEPTION(valuesAndValidatorIndex < 0,
    MissingValuesAndValidatorsTagException,
    "Error: All StringValidatorDependencies must have a " << 
    getValuesAndValidatorsTag() << "tag!" << std::endl << std::endl);

  XMLObject valuesAndValidatorTag = xmlObj.getChild(valuesAndValidatorIndex);
  for(int i=0; i < valuesAndValidatorTag.numChildren(); ++i){
    XMLObject child = valuesAndValidatorTag.getChild(i);
    std::string value = child.getRequired(getValueAttributeName());
    ParameterEntryValidator::ValidatorID valiID = 
      child.getRequired<ParameterEntryValidator::ValidatorID>(
        getValidatorIdAttributeName());
    TEST_FOR_EXCEPTION(validatorIDsMap.find(valiID) == validatorIDsMap.end(),
      MissingValidatorException,
      "Could not find a validator corresponding to the ID " << valiID <<
      " in the given validatorIDsMap!" << std::endl << std::endl);
    RCP<ParameterEntryValidator> validator = 
      validatorIDsMap.find(valiID)->second;
    valueValidatorMap.insert(
      StringValidatorDependency::ValueToValidatorPair(value, validator));
  }

  RCP<ParameterEntryValidator> defaultValidator = null;
  if(xmlObj.hasAttribute(getDefaultValidatorIDAttributeName())){
    ParameterEntryValidator::ValidatorID defaultValiID = 
      xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
        getDefaultValidatorIDAttributeName());
    TEST_FOR_EXCEPTION(
      validatorIDsMap.find(defaultValiID) == validatorIDsMap.end(),
      MissingValidatorException,
      "Could not find a validator (for the default validator) " <<
      "corresponding to the ID " << defaultValiID << 
      " in the given validatorIDsMap!" << std::endl << std::endl);
    defaultValidator = validatorIDsMap.find(defaultValiID)->second;
  }

  return rcp(new StringValidatorDependency(
    dependee, dependents, valueValidatorMap, defaultValidator));
}

void
BoolValidatorDependencyXMLConverter::convertSpecialValidatorAttributes(
  RCP<const ValidatorDependency> dependency,
  XMLObject& xmlObj,
  ValidatortoIDMap& validatorIDsMap) const
{
  RCP<const BoolValidatorDependency> castedDependency = 
    rcp_dynamic_cast<const BoolValidatorDependency>(dependency, true);  

  RCP<const ParameterEntryValidator> trueVali = 
    castedDependency->getTrueValidator();
  RCP<const ParameterEntryValidator> falseVali = 
    castedDependency->getFalseValidator();

  if(nonnull(trueVali)){
    if(validatorIDsMap.find(castedDependency->getTrueValidator()) == 
      validatorIDsMap.end()){
      validatorIDsMap.insert(castedDependency->getTrueValidator());
    }
    xmlObj.addAttribute(
      getTrueValidatorIDAttributeName(),
      validatorIDsMap.find(castedDependency->getTrueValidator())->second);
  }

  if(nonnull(falseVali)){
    if(validatorIDsMap.find(falseVali) == 
      validatorIDsMap.end()){
      validatorIDsMap.insert(falseVali);
    }
    xmlObj.addAttribute(
      getFalseValidatorIDAttributeName(),
      validatorIDsMap.find(falseVali)->second);
  }

}

RCP<ValidatorDependency> 
BoolValidatorDependencyXMLConverter::convertSpecialValidatorAttributes(
  const XMLObject& xmlObj,
  RCP<const ParameterEntry> dependee,
  const Dependency::ParameterEntryList dependents,
  const IDtoValidatorMap& validatorIDsMap) const
{

  RCP<ParameterEntryValidator> trueValidator = null;
  RCP<ParameterEntryValidator> falseValidator = null;

  if(xmlObj.hasAttribute(getTrueValidatorIDAttributeName())){

    ParameterEntryValidator::ValidatorID trueID = 
      xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
        getTrueValidatorIDAttributeName());
  
    TEST_FOR_EXCEPTION(
      validatorIDsMap.find(trueID)
      == 
      validatorIDsMap.end(),
      MissingValidatorException,
      "Could not find a Validator for the True validator " <<
      "with ID " << trueID <<
      " in the given validatorIDsMap!" << std::endl << std::endl);
  
    trueValidator = 
      validatorIDsMap.find(trueID)->second;
  }


  if(xmlObj.hasAttribute(getFalseValidatorIDAttributeName())){
    ParameterEntryValidator::ValidatorID falseID = 
      xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
        getFalseValidatorIDAttributeName());
  
    TEST_FOR_EXCEPTION(
      validatorIDsMap.find(falseID)
      == 
      validatorIDsMap.end(),
      MissingValidatorException,
      "Could not find a Validator for the False validator " <<
      "with ID " << falseID <<
      " in the given validatorIDsMap!" << std::endl << std::endl);
  
    falseValidator = 
      validatorIDsMap.find(falseID)->second;
  }

  return rcp(new BoolValidatorDependency(
    dependee, dependents, trueValidator, falseValidator));
}



} //namespace Teuchos

