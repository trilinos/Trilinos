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



#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#ifdef TEUCHOS_HAVE_QT
  #include <QDir>
#endif


namespace Teuchos{


VisualDependency::VisualDependency(
  RCP<const ParameterEntry> dependee, 
  RCP<ParameterEntry> dependent,
  bool showIf):
  Dependency(dependee, dependent),
  showIf_(showIf){}

VisualDependency::VisualDependency(
  RCP<const ParameterEntry> dependee,
  ParameterEntryList dependents,
  bool showIf):
  Dependency(dependee, dependents),
  showIf_(showIf){}

VisualDependency::VisualDependency(
  ConstParameterEntryList dependees, 
  RCP<ParameterEntry> dependent,
  bool showIf):
  Dependency(dependees, dependent),
  showIf_(showIf){}

VisualDependency::VisualDependency(
  ConstParameterEntryList dependees,
  ParameterEntryList dependents,
  bool showIf):
  Dependency(dependees, dependents),
  showIf_(showIf){}

bool VisualDependency::isDependentVisible() const{
  return dependentVisible_;
}

bool VisualDependency::getShowIf() const{
  return showIf_;
}

void VisualDependency::evaluate(){
  if((getDependeeState() && showIf_) || (!getDependeeState() && !showIf_)){
    dependentVisible_ = true;
  }
  else{
    dependentVisible_ = false;
  }
}
  
ValidatorDependency::ValidatorDependency( 
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent):
  Dependency(dependee, dependent){}

ValidatorDependency::ValidatorDependency(
  RCP<const ParameterEntry> dependee, 
  ParameterEntryList dependents):
  Dependency(dependee, dependents){}

StringVisualDependency::StringVisualDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  std::string value,
  bool showIf):
  VisualDependency(dependee, dependent, showIf),
  values_(ValueList(1,value))
{
  validateDep();
}

StringVisualDependency::StringVisualDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  const ValueList& values,
  bool showIf):
  VisualDependency(dependee, dependent, showIf),
  values_(values)
{
  validateDep();
}

StringVisualDependency::StringVisualDependency(
  RCP<const ParameterEntry> dependee, 
  Dependency::ParameterEntryList dependents, 
  const std::string& value,
  bool showIf):
  VisualDependency(dependee, dependents, showIf),
  values_(ValueList(1,value))
{
  validateDep();
}

StringVisualDependency::StringVisualDependency(
  RCP<const ParameterEntry> dependee, 
  Dependency::ParameterEntryList dependents, 
  const ValueList& values,
  bool showIf):
  VisualDependency(dependee, dependents, showIf),
  values_(values)
{
  validateDep();
}

const StringVisualDependency::ValueList& 
  StringVisualDependency::getValues() const
{
  return values_;
}

bool StringVisualDependency::getDependeeState() const{
  return find(values_.begin(), values_.end(), 
    getFirstDependeeValue<std::string>()) != values_.end();
}

std::string StringVisualDependency::getTypeAttributeValue() const{
  return "StringVisualDependency";
}

void StringVisualDependency::validateDep() const{
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<std::string>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "String Visual Dependency must be of type " 
    << TypeNameTraits<std::string>::name() << std::endl <<
    "Type encountered: " << getFirstDependee()->getAny().typeName() << 
    std::endl << std::endl);
}

RCP<StringVisualDependency> 
  DummyObjectGetter<StringVisualDependency>::getDummyObject()
{
  std::string blahString = "blah";
  return rcp(new StringVisualDependency(
    rcp(new ParameterEntry(blahString)),
    DummyObjectGetter<ParameterEntry>::getDummyObject(),
    "i'm a dummy"));
}

BoolVisualDependency::BoolVisualDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  bool showIf):
  VisualDependency(dependee, dependent, showIf)
{
  validateDep();
}

BoolVisualDependency::BoolVisualDependency(
  RCP<const ParameterEntry> dependee, 
  Dependency::ParameterEntryList dependents, 
  bool showIf):
  VisualDependency(dependee, dependents, showIf)
{
  validateDep();
}

bool BoolVisualDependency::getDependeeState() const{
  return getFirstDependeeValue<bool>();
}

std::string BoolVisualDependency::getTypeAttributeValue() const{
  return "BoolVisualDependency";
}

void BoolVisualDependency::validateDep() const{
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<bool>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "Bool Visual Dependency must be of type " << 
    TypeNameTraits<bool>::name() << std::endl <<
    "Encountered type: " << getFirstDependee()->getAny().typeName() << 
    std::endl << std::endl);
}

RCP<BoolVisualDependency> 
  DummyObjectGetter<BoolVisualDependency>::getDummyObject()
{
  return rcp(new BoolVisualDependency(
    rcp(new ParameterEntry(true)),
    DummyObjectGetter<ParameterEntry>::getDummyObject()));
}

ConditionVisualDependency::ConditionVisualDependency(
  RCP<const Condition> condition,
  RCP<ParameterEntry> dependent, 
  bool showIf):
  VisualDependency(condition->getAllParameters(), dependent, showIf),
  condition_(condition)
{
  validateDep();
}

ConditionVisualDependency::ConditionVisualDependency(
  RCP<const Condition> condition,
  Dependency::ParameterEntryList dependents,
  bool showIf):
  VisualDependency(condition->getAllParameters(), dependents, showIf),
  condition_(condition)
{
  validateDep();
}

RCP<const Condition> ConditionVisualDependency::getCondition() const{
  return condition_;
}

bool ConditionVisualDependency::getDependeeState() const{
  return condition_->isConditionTrue();
}

std::string ConditionVisualDependency::getTypeAttributeValue() const{
  return "ConditionVisualDependency";
}

RCP<ConditionVisualDependency> 
  DummyObjectGetter<ConditionVisualDependency>::getDummyObject()
{
  return rcp(new ConditionVisualDependency(
    DummyObjectGetter<NotCondition>::getDummyObject(),
    DummyObjectGetter<ParameterEntry>::getDummyObject()));
}

StringValidatorDependency::StringValidatorDependency(
  RCP<const ParameterEntry> dependee, 
  RCP<ParameterEntry> dependent,
  ValueToValidatorMap valuesAndValidators, 
  RCP<ParameterEntryValidator> defaultValidator):
  ValidatorDependency(dependee, dependent),
  valuesAndValidators_(valuesAndValidators),
  defaultValidator_(defaultValidator)
{
  validateDep();
}

StringValidatorDependency::StringValidatorDependency(
  RCP<const ParameterEntry> dependee, 
  Dependency::ParameterEntryList dependents,
  ValueToValidatorMap valuesAndValidators, 
  RCP<ParameterEntryValidator> defaultValidator):
  ValidatorDependency(dependee, dependents),
  valuesAndValidators_(valuesAndValidators),
  defaultValidator_(defaultValidator)
{
  validateDep();
}

const StringValidatorDependency::ValueToValidatorMap& 
  StringValidatorDependency::getValuesAndValidators() const
{
  return valuesAndValidators_;
}

RCP<const ParameterEntryValidator> 
  StringValidatorDependency::getDefaultValidator() const
{
  return defaultValidator_;
}

void StringValidatorDependency::evaluate(){
  std::string currentDependeeValue = getFirstDependeeValue<std::string>();
  for(
    ParameterEntryList::iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    if(
      valuesAndValidators_.find(currentDependeeValue) 
      == 
      valuesAndValidators_.end())
    {
      (*it)->setValidator(defaultValidator_);
    }
    else{
      (*it)->setValidator(valuesAndValidators_[currentDependeeValue]);
    }
  }
}

std::string StringValidatorDependency::getTypeAttributeValue() const{
  return "StringValidatorDependency";
}

void StringValidatorDependency::validateDep() const{
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<std::string>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "String Validator Dependency must be of type " <<
    TypeNameTraits<std::string>::name() << std::endl <<
    "Type Encountered: " << getFirstDependee()->getAny().typeName() <<
    std::endl << std::endl);

  TEST_FOR_EXCEPTION(
    valuesAndValidators_.size() < 1,
    InvalidDependencyException,
    "The valuesAndValidatord map for a string validator dependency must "
    "have at least one entry!" << std::endl << std::endl);
  ValueToValidatorMap::const_iterator it = valuesAndValidators_.begin();
  RCP<const ParameterEntryValidator> firstVali = (it->second);
  ++it;
  for(; it != valuesAndValidators_.end(); ++it){
    TEST_FOR_EXCEPTION( typeid(*firstVali) != typeid(*(it->second)),
      InvalidDependencyException,
      "Ay no! All of the validators in a StringValidatorDependency "
      "must have the same type.");
   }
}


RCP<StringValidatorDependency >
  DummyObjectGetter<StringValidatorDependency>::getDummyObject()
{
  std::string blahString = "blah";
  StringValidatorDependency::ValueToValidatorMap dummyMap;
  dummyMap.insert(StringValidatorDependency::ValueToValidatorPair("blah",
    DummyObjectGetter<FileNameValidator>::getDummyObject()));
  return rcp(new StringValidatorDependency(
    rcp(new ParameterEntry(blahString)),
    DummyObjectGetter<ParameterEntry>::getDummyObject(),
    dummyMap));
}
  
BoolValidatorDependency::BoolValidatorDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  RCP<const ParameterEntryValidator> trueValidator,
  RCP<const ParameterEntryValidator> falseValidator):
  ValidatorDependency(dependee, dependent),
  trueValidator_(trueValidator),
  falseValidator_(falseValidator)
{
  validateDep();
}

BoolValidatorDependency::BoolValidatorDependency(
  RCP<const ParameterEntry> dependee,
  Dependency::ParameterEntryList dependents,
  RCP<const ParameterEntryValidator> trueValidator,
  RCP<const ParameterEntryValidator> falseValidator):
  ValidatorDependency(dependee, dependents),
  trueValidator_(trueValidator),
  falseValidator_(falseValidator)
{
  validateDep();
}

void BoolValidatorDependency::evaluate(){
  bool dependeeValue = getFirstDependeeValue<bool>();
  for(
    ParameterEntryList::iterator it = getDependents().begin();
    it != getDependents().end(); 
    ++it)
  { 
    dependeeValue ? 
      (*it)->setValidator(trueValidator_) 
      :
      (*it)->setValidator(falseValidator_);
  }
}

RCP<const ParameterEntryValidator> 
  BoolValidatorDependency::getTrueValidator() const
{
  return trueValidator_;
}

RCP<const ParameterEntryValidator> 
  BoolValidatorDependency::getFalseValidator() const
{
  return falseValidator_;
}

std::string BoolValidatorDependency::getTypeAttributeValue() const{
  return "BoolValidatorDependency";
}
  

void BoolValidatorDependency::validateDep() const{

  #ifdef TEUCHOS_DEBUG
  if(falseValidator_.is_null() && trueValidator_.is_null()){
    std::cerr << "WARNING: Both the true and false validators " <<
    "for a BoolValidatorDependency seem to be null. You sure you want that?" <<
    std::endl << std::endl;
  }
  #endif

  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<bool>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "Bool Validator Dependency must be of type " <<
    TypeNameTraits<bool>::name() << std::endl <<
    "Encountered type: " << getFirstDependee()->getAny().typeName() <<
    std::endl << std::endl);

  if(!falseValidator_.is_null() && !trueValidator_.is_null()){
    TEST_FOR_EXCEPTION(typeid(*falseValidator_) != typeid(*trueValidator_),
      InvalidDependencyException,
      "Ay no! The true and false validators of a Bool Validator Dependency "
      "must be the same type! " <<std::endl << std::endl);
  }
  
}

RCP<BoolValidatorDependency>
  DummyObjectGetter<BoolValidatorDependency>::getDummyObject()
{
  return rcp(new BoolValidatorDependency(
    rcp(new ParameterEntry(true)),
    DummyObjectGetter<ParameterEntry>::getDummyObject(),
    null, null));
}

}

