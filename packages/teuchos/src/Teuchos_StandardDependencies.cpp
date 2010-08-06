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

VisualDependency::VisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf)
:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList), showIf_(showIf){}

VisualDependency::VisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
ParameterParentMap dependents, bool showIf)
:Dependency(dependeeName, dependeeParentList, dependents), showIf_(showIf){}

VisualDependency::VisualDependency(ParameterParentMap dependees, std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf)
:Dependency(dependees, dependentName, dependentParentList), showIf_(showIf){}

VisualDependency::VisualDependency(ParameterParentMap dependees, ParameterParentMap dependents, bool showIf) 
:Dependency(dependees, dependents), showIf_(showIf){}

bool VisualDependency::isDependentVisible() const{
  return dependentVisible_;
}

ValidatorDependency::ValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList)
:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList){}

ValidatorDependency::ValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
ParameterParentMap dependents)
:Dependency(dependeeName, dependeeParentList, dependents){}


StringVisualDependency::StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList, std::string value, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList, showIf), values_(ValueList(1,value)){
  validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, std::string dependentName, 
RCP<ParameterList> parentList, std::string value, bool showIf)
:VisualDependency(dependeeName, parentList, dependentName, parentList, showIf), values_(ValueList(1,value)){
  validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList, const ValueList& values, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList, showIf), values_(values){
  validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, std::string dependentName, 
RCP<ParameterList> parentList, const ValueList& values, bool showIf)
:VisualDependency(dependeeName, parentList, dependentName, parentList, showIf), values_(values){
  validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList, 
ParameterParentMap dependents, std::string value, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependents, showIf), values_(ValueList(1,value)){
  validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
ParameterParentMap dependents, const ValueList& values, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependents, showIf), values_(values){
  validateDep();
}

bool StringVisualDependency::getDependeeState() const{
  return find(values_.begin(), values_.end(), getFirstDependeeValue<std::string>()) != values_.end();
}

void StringVisualDependency::validateDep() const{
  /*
   * This error should never get thrown, unless someone
   * is doing something wonky in a sublcass.
   */
  TEST_FOR_EXCEPTION(getDependees().size() != 1,
    InvalidDependencyException,
    "Uh oh. Looks like you tried to make a "
    "String Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
    "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
    "be doing some subclassing you sly-dog ;)\n\n" 
    "Error: A String Visual Dependency must have exactly 1 dependee. " 
    "You have tried to assign it "<< getDependees().size() << " dependees.\n" 
    "Dependees: " << getDependeeNamesString() << "\n" 
    "Dependents: " << getDependentNamesString());

  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<std::string>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "String Visual Dependency must be of type string.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependent: " << getDependentNamesString());

}

BoolVisualDependency::BoolVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList, showIf){
  validateDep();
}

BoolVisualDependency::BoolVisualDependency(std::string dependeeName, std::string dependentName, 
RCP<ParameterList> parentList, bool showIf)
:VisualDependency(dependeeName, parentList, dependentName, parentList, showIf){
  validateDep();
}

BoolVisualDependency::BoolVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
ParameterParentMap dependents, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependents, showIf){
  validateDep();
}

bool BoolVisualDependency::getDependeeState() const{
  return getFirstDependeeValue<bool>();
}

void BoolVisualDependency::validateDep() const{
  /*
   * This error should never get thrown, unless someone
   * is doing something wonky in a sublcass.
   */
  TEST_FOR_EXCEPTION(getDependees().size() != 1,
    InvalidDependencyException,
    "Uh oh. Looks like you tried to make a "
    "Bool Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
    "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
    "be doing some subclassing you sly-dog ;)\n\n" 
    "Error: A Bool Visual Dependency must have exactly 1 dependee. " 
    "You have tried to assign it "<< getDependees().size() << " dependees.\n" 
    "Dependees: " << getDependeeNamesString() << "\n" 
    "Dependents: " << getDependentNamesString());

  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<bool>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "Bool Visual Dependency must be of type bool.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependents: " << getDependentNamesString());

}

ConditionVisualDependency::ConditionVisualDependency(RCP<Condition> condition,
std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf)
:VisualDependency(condition->getAllParameters(), dependentName, dependentParentList, showIf), condition_(condition){
  validateDep();
}

ConditionVisualDependency::ConditionVisualDependency(RCP<Condition> condition, ParameterParentMap dependents, bool showIf)
:VisualDependency(condition->getAllParameters(), dependents, showIf), condition_(condition){
  validateDep();
}

bool ConditionVisualDependency::getDependeeState() const{
  return condition_->isConditionTrue();
}

void ConditionVisualDependency::validateDep() const {}

NumberArrayLengthDependency::NumberArrayLengthDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList, int (*func)(int))
:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList), func_(func){
  validateDep();
}

NumberArrayLengthDependency::NumberArrayLengthDependency(std::string dependeeName, std::string dependentName, 
RCP<ParameterList> parentList, int (*func)(int))
:Dependency(dependeeName, parentList, dependentName, parentList), func_(func){
  validateDep();
}

NumberArrayLengthDependency::NumberArrayLengthDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
ParameterParentMap dependents, int (*func)(int))
:Dependency(dependeeName, dependeeParentList, dependents), func_(func){
  validateDep();
}

int NumberArrayLengthDependency::runFunction(int argument) const{
  if(func_ !=0){
    return (*func_)(argument);
  }
  else{
    return argument;
  }
}

template <class T>
void NumberArrayLengthDependency::modifyArrayLength(int newLength, ParameterEntry* dependentToModify){
  const Array<T> originalArray = any_cast<Array<T> >(dependentToModify->getAny()); 
  RCP<const EnhancedNumberValidator<T> > potentialValidator = null;
  if(!is_null(dependentToModify->validator())){
      potentialValidator = rcp_dynamic_cast<const ArrayNumberValidator<T> >(dependentToModify->validator(),true)->getPrototype();
  }
  Array<T> newArray;
  int i;
  for(i=0; i<originalArray.size() && i<newLength; ++i){
    newArray.append(originalArray[i]);
  }
  for(;i<newLength;++i){
    if(is_null(potentialValidator)){
      newArray.append(0);
    }
    else{
      newArray.append(potentialValidator->getMin());
    }
  }
  dependentToModify->setValue(newArray, false, dependentToModify->docString(), dependentToModify->validator());
}

template<>
void NumberArrayLengthDependency::modifyArrayLength<std::string>(int newLength, ParameterEntry* dependentToModify){
  const Array<std::string> originalArray = any_cast<Array<std::string> >(dependentToModify->getAny()); 
  Array<std::string> newArray;
  RCP<const ParameterEntryValidator> validator = dependentToModify->validator();
  int i;
  for(i=0; i<originalArray.size() && i<newLength; ++i){
    newArray.append(originalArray[i]);
  }
  for(;i<newLength;++i){
    if(is_null(validator)){
      newArray.append(" ");
    }
    else if(!is_null(rcp_dynamic_cast<const ArrayFileNameValidator>(validator))){
      #ifdef TEUCHOS_HAVE_QT
        newArray.append(QDir::homePath().toStdString());
      #else
        newArray.append("");
      #endif
    }
    else{
      newArray.append(validator->validStringValues()->at(0));
    }
  }
  dependentToModify->setValue(newArray, false, dependentToModify->docString(), dependentToModify->validator());
}

void NumberArrayLengthDependency::evaluate(){
  int newLength = runFunction(getFirstDependeeValue<int>());
  TEST_FOR_EXCEPTION(newLength <0, Exceptions::InvalidParameterValue,
    "Ruh Roh Shaggy! Looks like a dependency tried to set the length of the Array(s) in the " <<
    getDependentNamesString() << " parameter(s) to a negative number. Silly. You can't have an Array " <<
    "with a negative length! You should probably contact the maintainer of this program and " <<
    "give him or her the following information: \n\n" <<
    "Error:\n" <<
    "An attempt was made to set the length of an Array to a negative number by a NumberArrayLengthDependency\n" <<
    "Dependency Type: " << getType() << "\n" <<
    "Problem Dependee: " << getFirstDependeeName() <<
    "Problem Dependents: " << getDependentNamesString());
  ParameterEntry *currentDependent;
  for(ParameterParentMap::const_iterator it = getDependents().begin(); it != getDependents().end(); ++it){ 
    currentDependent = it->second->getEntryPtr(it->first);
    if(currentDependent->getAny().type() == typeid(Array<int>)){
      modifyArrayLength<int>(newLength, currentDependent);
    }
    else if(currentDependent->getAny().type() == typeid(Array<short>)){
      modifyArrayLength<short>(newLength, currentDependent);
    }
    else if(currentDependent->getAny().type() == typeid(Array<double>)){
      modifyArrayLength<double>(newLength, currentDependent);
    }
    else if(currentDependent->getAny().type() == typeid(Array<float>)){
      modifyArrayLength<float>(newLength, currentDependent);
    }
    else if(currentDependent->getAny().type() == typeid(Array<std::string>)){
      modifyArrayLength<std::string>(newLength, currentDependent);
    }
  }
}

void NumberArrayLengthDependency::validateDep() const{
  /*
   * This error should never get thrown, unless someone
   * is doing something wonky in a sublcass.
   */
  TEST_FOR_EXCEPTION(getDependees().size() != 1,
    InvalidDependencyException,
    "Uh oh. Looks like you tried to make a "
    "Number Array Length Dependency doesn't have exactly one currentDependee. This is kind of a problem. " 
    "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
    "be doing some subclassing you sly-dog ;)\n\n" 
    "Error: A Number Array Length Dependency must have exactly 1 currentDependee. " 
    "You have tried to assign it "<< getDependees().size() << " dependees.\n" 
    "Dependees: " << getDependeeNamesString() << "\n" 
    "Dependents: " << getDependentNamesString());

  #ifndef HAVE_TEUCHOS_LONG_LONG_INT
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<int>() 
    && !getFirstDependee()->isType<short>() 
    && !getFirstDependee()->isType<long>() ,
    InvalidDependencyException,
    "Ay no! The currentDependee for an "
    "Array Length Dependency must be either of type short or of type int.\n" 
    "Problem currentDependee: " << getFirstDependeeName() << "\n" 
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependent: " << getDependentNamesString());
  #else
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<int>() 
    && !getFirstDependee()->isType<short>() 
    && !getFirstDependee()->isType<long>() 
    && !getFirstDependee()->isType<long long>(),
    InvalidDependencyException,
    "Ay no! The currentDependee for an "
    "Array Length Dependency must be either of type short or of type int.\n" 
    "Problem currentDependee: " << getFirstDependeeName() << "\n" 
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependent: " << getDependentNamesString());
  #endif
  ParameterEntry *currentDependent;
  for(ParameterParentMap::const_iterator it = getDependents().begin(); it != getDependents().end(); ++it){ 
    currentDependent = it->second->getEntryPtr(it->first);
    TEST_FOR_EXCEPTION(!currentDependent->isArray(),
      InvalidDependencyException,
      "Ay no! The dependent of an "
      "Array Length Dependency must be an array.\n"
      "Problem dependent: " << it->first << "\n" 
      "Actual type: " << currentDependent->getAny().typeName() << "\n"
      "Dependees: " << getDependeeNamesString());
  }
}

StringValidatorDependency::StringValidatorDependency(std::string currentDependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList,
ValueToValidatorMap valuesAndValidators, RCP<ParameterEntryValidator> defaultValidator)
:ValidatorDependency(currentDependeeName, dependeeParentList, dependentName, dependentParentList), defaultValidator_(defaultValidator), valuesAndValidators_(valuesAndValidators){
  validateDep();
}

StringValidatorDependency::StringValidatorDependency(std::string currentDependeeName, std::string dependentName, 
RCP<ParameterList> parentList, ValueToValidatorMap valuesAndValidators, 
RCP<ParameterEntryValidator> defaultValidator)
:ValidatorDependency(currentDependeeName, parentList, dependentName, parentList), defaultValidator_(defaultValidator), valuesAndValidators_(valuesAndValidators){
  validateDep();
}

StringValidatorDependency::StringValidatorDependency(std::string currentDependeeName, RCP<ParameterList> dependeeParentList,
ParameterParentMap dependents,
ValueToValidatorMap valuesAndValidators, RCP<ParameterEntryValidator> defaultValidator)
:ValidatorDependency(currentDependeeName, dependeeParentList, dependents), defaultValidator_(defaultValidator), valuesAndValidators_(valuesAndValidators){
  validateDep();
}

void StringValidatorDependency::evaluate(){
  std::string currentDependeeValue = getFirstDependeeValue<std::string>();
  ParameterEntry *currentDependent;
  for(ParameterParentMap::const_iterator it = getDependents().begin(); it != getDependents().end(); ++it){ 
    currentDependent = it->second->getEntryPtr(it->first);
    if(valuesAndValidators_.find(currentDependeeValue) == valuesAndValidators_.end()){
      currentDependent->setValidator(defaultValidator_);
    }
    else{
      currentDependent->setValidator(valuesAndValidators_[currentDependeeValue]);
    }
  }
}
void StringValidatorDependency::validateDep() const{
  /*
   * This error should never get thrown, unless someone
   * is doing something wonky in a sublcass.
   */
  TEST_FOR_EXCEPTION(getDependees().size() != 1,
    InvalidDependencyException,
    "Uh oh. Looks like you tried to make a "
    "String Validator Dependency doesn't have exactly one dependee. This is kind of a problem. " 
    "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
    "be doing some subclassing you sly-dog ;)\n\n" 
    "Error: A String Validator Dependency must have exactly 1 dependee. "
    "You have tried to assign it "<< getDependees().size() << " dependees.\n"
    "Dependees: " << getDependeeNamesString() << "\n"
    "Dependents: " << getDependentNamesString());

  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<std::string>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "String Validator Dependency must be of type string.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependent: " << getDependentNamesString());

  for(ValueToValidatorMap::const_iterator it = valuesAndValidators_.begin(); it != valuesAndValidators_.end(); ++it){
    ParameterEntry *currentDependent;
    for(ParameterParentMap::const_iterator it2 = getDependents().begin(); it2 != getDependents().end(); ++it2){ 
      currentDependent = it2->second->getEntryPtr(it2->first);
      TEST_FOR_EXCEPTION(typeid(*(currentDependent->validator().get())) != typeid(*(it->second.get())),
        InvalidDependencyException,
        "Ay no! The validator of a dependent of a "
        "String Validator Dependency must be the same type as all of the validators "
        "in the valuesAndValidators map. "
        "Note this means that the dependent must have an initial validator.\n"
        "Problem dependent: " << it->first << "\n"
        "Validator Type: " << typeid(*(currentDependent->validator())).name() << "\n"
        "Problem Map Value: " << it->first << "\n"
        "The Validator type associated with the Map Value: " << typeid(*(it->second)).name());
    }
  }
}

BoolValidatorDependency::BoolValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
std::string dependentName, RCP<ParameterList> dependentParentList,
RCP<const ParameterEntryValidator> trueValidator, RCP<const ParameterEntryValidator> falseValidator)
:ValidatorDependency(dependeeName, dependeeParentList, dependentName, dependentParentList), trueValidator_(trueValidator), falseValidator_(falseValidator){
  validateDep();
}

BoolValidatorDependency::BoolValidatorDependency(std::string dependeeName, std::string dependentName, 
RCP<ParameterList> parentList, RCP<const ParameterEntryValidator> trueValidator,
RCP<const ParameterEntryValidator> falseValidator)
:ValidatorDependency(dependeeName, parentList, dependentName, parentList), trueValidator_(trueValidator), falseValidator_(falseValidator){
  validateDep();
}

BoolValidatorDependency::BoolValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
ParameterParentMap dependents,
RCP<const ParameterEntryValidator> trueValidator, RCP<const ParameterEntryValidator> falseValidator)
:ValidatorDependency(dependeeName, dependeeParentList, dependents), trueValidator_(trueValidator), falseValidator_(falseValidator){
  validateDep();
}

void BoolValidatorDependency::evaluate(){
  bool dependeeValue = getFirstDependeeValue<bool>();
  ParameterEntry *currentDependent;
  for(ParameterParentMap::const_iterator it = getDependents().begin(); it != getDependents().end(); ++it){ 
    currentDependent = it->second->getEntryPtr(it->first);
    dependeeValue ? currentDependent->setValidator(trueValidator_) : currentDependent->setValidator(falseValidator_);
  }
}

void BoolValidatorDependency::validateDep() const{
  /*
   * This error should never get thrown, unless someone
   * is doing something wonky in a sublcass.
   */
  TEST_FOR_EXCEPTION(getDependees().size() != 1,
    InvalidDependencyException,
    "Uh oh. Looks like you tried to make a "
    "Bool Validator Dependency doesn't have exactly one dependee. This is kind of a problem. " 
    "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
    "be doing some subclassing you sly-dog ;)\n\n" 
    "Error: A Bool Validator Dependency must have exactly 1 dependee. "
    "You have tried to assign it "<< getDependees().size() << " dependees.\n"
    "Dependees: " << getDependeeNamesString() << "\n" 
    "Dependents: " << getDependentNamesString());

  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<bool>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "Bool Validator Dependency must be of type boolean.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependent: " << getDependentNamesString());
  
  ParameterEntry *currentDependent;
  for(ParameterParentMap::const_iterator it = getDependents().begin(); it != getDependents().end(); ++it){ 
    currentDependent = it->second->getEntryPtr(it->first);
    TEST_FOR_EXCEPTION(typeid(*(currentDependent->validator().get())) != typeid(*(trueValidator_.get())),
      InvalidDependencyException,
      "Ay no! The validator of a dependent of a "
      "Bool Validator Dependency must be the same type as the \"true\" validator. "
      "Note this means that the dependent must have an initial validator.\n"
      "Problem dependent: " << it->first << "\n"
      "Validator Type: " << typeid(*(currentDependent->validator().get())).name() << "\n"
      "Type of the \"true\" validator: " << typeid(*(trueValidator_.get())).name());

    TEST_FOR_EXCEPTION(typeid(*(currentDependent->validator().get())) != typeid(*(falseValidator_.get())),
      InvalidDependencyException,
      "Ay no! The validator of a dependent of a "
      "Bool Validator Dependency must be the same type as the \"false\" validator. "
      "Note this means that the dependent must have an initial validator.\n"
      "Problem dependent: " << it->first << "\n"
      "Validator Type: " << typeid(*(currentDependent->validator().get())).name() << "\n"
      "Type of the \"false\" validator: " << typeid(*(falseValidator_.get())).name());
  }
}

}

