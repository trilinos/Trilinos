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
  RCP<ParameterEntry> dependee, 
  RCP<ParameterEntry> dependent,
  bool showIf):
  Dependency(dependee, dependent),
  showIf_(showIf){}

VisualDependency::VisualDependency(
  RCP<ParameterEntry> dependee,
  DependentList& dependents,
  bool showIf):
  Dependency(dependee, dependents),
  showIf_(showIf){}

VisualDependency::VisualDependency(
  DependeeList& dependees, 
  RCP<ParameterEntry> dependent,
  bool showIf):
  Dependency(dependees, dependent),
  showIf_(showIf){}

VisualDependency::VisualDependency(
  DependeeList dependees,
  DependentList dependents,
  bool showIf):
  Dependency(dependees, dependents),
  showIf_(showIf){}



bool VisualDependency::isDependentVisible() const{
  return dependentVisible_;
}

void VisualDependency::evaluate(){
  if((getDependeeState() && showIf_) || (!getDependeeState() && !showIf_)){
    dependentVisible_ = true;
  }
  else{
    dependentVisible_ = false;
  }
}
  
ValidatorDependency::ValidatorDependency( RCP<ParameterEntry> dependee,
  RCP<ParameterEntry> dependent):
  Dependency(dependee, dependent){}

ValidatorDependency::ValidatorDependency(RCP<ParameterEntry> dependee, 
  DependentList dependents):
  Dependency(dependee, dependents){}

StringVisualDependency::StringVisualDependency(
  RCP<ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  std::string value,
  bool showIf):
  VisualDependency(dependee, dependent, showIf),
  values_(ValueList(1,value))
{
  validateDep();
}

StringVisualDependency::StringVisualDependency(
  RCP<ParameterEntry> dependent,
  RCP<ParameterEntry> dependee,
  const ValueList& values,
  bool showIf):
  VisualDependency(dependee, dependent, showIf),
  values_(values)
{
  validateDep();
}

StringVisualDependency::StringVisualDependency(
  RCP<ParameterEntry> dependee, 
  Dependency::DependentList& dependents, 
  const std::string& value,
  bool showIf):
  VisualDependency(dependee, dependents, showIf),
  values_(ValueList(1,value))
{
  validateDep();
}

StringVisualDependency::StringVisualDependency(
  RCP<ParameterEntry> dependee, 
  Dependency::DependentList& dependents, 
  const ValueList& values,
  bool showIf):
  VisualDependency(dependee, dependents, showIf),
  values_(values)
{
  validateDep();
}

void StringVisualDependency::validateDep() const{
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<std::string>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "String Visual Dependency must be of type string.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n");
}

BoolVisualDependency::BoolVisualDependency(
  RCP<ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  bool showIf):
  VisualDependency(dependee, dependent, showIf)
{
  validateDep();
}

BoolVisualDependency::BoolVisualDependency(
  RCP<ParameterEntry> dependee, 
  Dependency::DependentList dependents, 
  bool showIf):
  VisualDependency(dependee, dependent, showIf)
{
  validateDep();
}

void BoolVisualDependency::validateDep() const{
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<bool>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "Bool Visual Dependency must be of type bool.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependents: " << getDependentNamesString());
}

ConditionVisualDependency::ConditionVisualDependency(
  RCP<Condition> condition,
  RCP<ParameterEntry> dependent, 
  bool showIf):
  VisualDependency(condition->getAllParameters(), dependent, showIf),
  condition_(condition)
{
  validateDep();
}

ConditionVisualDependency::ConditionVisualDependency(
  RCP<Condition> condition,
  Dependency::DependentList dependents,
  bool showIf):
  VisualDependency(condition->getAllParameters(), dependents, showIf),
  condition_(condition)
{
  validateDep();
}

NumberArrayLengthDependency::NumberArrayLengthDependency(
  RCP<ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  int (*func)(int) = 0):
  Dependency(dependee, dependent),
  func_(func)
{
  validateDep();
}

NumberArrayLengthDependency::NumberArrayLengthDependency(
  RCP<ParameterEntry> dependee,
  DependentList& dependent,
  int (*func)(int) = 0):
  Dependency(dependee, dependents),
  func_(func)
{
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
void NumberArrayLengthDependency::modifyArrayLength(
  int newLength, ParameterEntry* dependentToModify)
{
  const Array<T> originalArray = 
    any_cast<Array<T> >(dependentToModify->getAny()); 
  RCP<const EnhancedNumberValidator<T> > potentialValidator = null;
  if(!is_null(dependentToModify->validator())){
      potentialValidator = rcp_dynamic_cast<const ArrayNumberValidator<T> >
        (dependentToModify->validator(),true)->getPrototype();
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
  dependentToModify->setValue(newArray,
    false, dependentToModify->docString(), dependentToModify->validator());
}

template<>
void NumberArrayLengthDependency::modifyArrayLength<std::string>(
  int newLength, ParameterEntry* dependentToModify)
{
  const Array<std::string> originalArray = 
    any_cast<Array<std::string> >(dependentToModify->getAny()); 
  Array<std::string> newArray;
  RCP<const ParameterEntryValidator> validator = 
    dependentToModify->validator();
  int i;
  for(i=0; i<originalArray.size() && i<newLength; ++i){
    newArray.append(originalArray[i]);
  }
  for(;i<newLength;++i){
    if(is_null(validator)){
      newArray.append(" ");
    }
    else if(
      !is_null(rcp_dynamic_cast<const ArrayFileNameValidator>(validator)))
    {
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
  dependentToModify->setValue(
    newArray, false, dependentToModify->docString(), 
    dependentToModify->validator());
}

void NumberArrayLengthDependency::evaluate(){
  int newLength = runFunction(getFirstDependeeValue<int>());
  TEST_FOR_EXCEPTION(newLength <0, Exceptions::InvalidParameterValue,
    "Ruh Roh Shaggy! Looks like a dependency tried to set the length "
    "of the Array(s) in the " << getDependentNamesString() << 
    " parameter(s) to a negative number. Silly. You can't have "
    "an Array with a negative length! You should probably contact the "
    "maintainer of this program and "
    "give him or her the following information: \n\n" <<
    "Error:\n" <<
    "An attempt was made to set the length of an Array to a negative "
    "number by a NumberArrayLengthDependency\n");
  ParameterEntry *currentDependent;
  for(
    ParameterParentMap::const_iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
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
    else if(
      currentDependent->getAny().type() == typeid(Array<std::string>))
    {
      modifyArrayLength<std::string>(newLength, currentDependent);
    }
  }
}

void NumberArrayLengthDependency::validateDep() const{
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<int>() 
    && !getFirstDependee()->isType<short>(),
    InvalidDependencyException,
    "Ay no! The Dependee for an "
    "Array Length Dependency must be either of type short or "
    "of type int.\n" 
    "Problem Dependee type: " << getFirstDependee()->getAny().typeName() 
    << "\n");

  ParameterEntry *currentDependent;
  for(
    ParameterParentMap::const_iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    currentDependent = it->second->getEntryPtr(it->first);
    TEST_FOR_EXCEPTION(!currentDependent->isArray(),
      InvalidDependencyException,
      "Ay no! The dependent of an "
      "Array Length Dependency must be an array.\n"
      "Actual Dependent type: "
      << currentDependent->getAny().typeName() << "\n");
  }
}

StringValidatorDependency(
  RCP<ParameterEntry> dependee, 
  RCP<ParameterEntry> dependent,
  ValueToValidatorMap valuesAndValidators, 
  RCP<ParameterEntryValidator> defaultValidator):
  ValidatorDependency(dependee, dependent),
  valuesAndValidators_(valuesAndValidators),
  defaultValidator_(defaultValidator)
{
  validateDep();
}

StringValidatorDependency(
  RCP<ParameterEntry> dependeeName, 
  Dependency::DependeeList& dependents,
  ValueToValidatorMap valuesAndValidators, 
  RCP<ParameterEntryValidator> defaultValidator):
  ValidatorDependency(dependee, dependent),
  valuesAndValidators_(valuesAndValidators),
  defaultValidator_(defaultValidator)
{
  validateDep();
}

void StringValidatorDependency::evaluate(){
  std::string currentDependeeValue = getFirstDependeeValue<std::string>();
  ParameterEntry *currentDependent;
  for(
    ParameterParentMap::const_iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    currentDependent = it->second->getEntryPtr(it->first);
    if(
      valuesAndValidators_.find(currentDependeeValue) 
      == 
      valuesAndValidators_.end())
    {
      currentDependent->setValidator(defaultValidator_);
    }
    else{
      currentDependent->setValidator(
        valuesAndValidators_[currentDependeeValue]);
    }
  }
}
void StringValidatorDependency::validateDep() const{
  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<std::string>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "String Validator Dependency must be of type string.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependent: " << getDependentNamesString());

  ValueToValidatorMap::const_iterator it = valuesAndValidators_.begin();
  type_info firstValidatorType = typeid(it->second);
  ++it;
  for(; it != valuesAndValidators_.end(); ++it){
    TEST_FOR_EXCEPTION( firstValidatorType != typeid(it->second),
      InvalidDependencyException,
      "Ay no! All of the validators in a StringValidatorDependency"
      "must have the same type.");
   }
}

BoolValidatorDependency(
  RCP<ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  RCP<const ParameterEntryValidator> trueValidator,
  RCP<const ParameterEntryValidator> falseValidator):
  ValidatorDependency(dependee, dependent),
  trueValidator_(trueValidator);
  falseValidator_(falseValidator);

BoolValidatorDependency(
  RCP<ParameterEntry> dependee,
  Dependency::DependentList dependents,
  RCP<const ParameterEntryValidator> trueValidator,
  RCP<const ParameterEntryValidator> falseValidator):
  ValidatorDependency(dependee, dependents),
  trueValidator_(trueValidator);
  falseValidator_(falseValidator);

void BoolValidatorDependency::evaluate(){
  bool dependeeValue = getFirstDependeeValue<bool>();
  ParameterEntry *currentDependent;
  for(
    ParameterParentMap::const_iterator it = getDependents().begin();
    it != getDependents().end(); 
    ++it)
  { 
    currentDependent = it->second->getEntryPtr(it->first);
    dependeeValue ? 
      currentDependent->setValidator(trueValidator_) 
      :
      currentDependent->setValidator(falseValidator_);
  }
}

void BoolValidatorDependency::validateDep() const{

  TEST_FOR_EXCEPTION(!getFirstDependee()->isType<bool>(),
    InvalidDependencyException,
    "Ay no! The dependee of a "
    "Bool Validator Dependency must be of type boolean.\n"
    "Problem dependee: " << getFirstDependeeName() << "\n"
    "Actual type: " << getFirstDependee()->getAny().typeName() << "\n"
    "Dependent: " << getDependentNamesString());

  TEST_FOR_EXCEPTION(typeid(falseValidator_) != typeid(trueValidator_),
    InvalidDependencyException,
    "Ay no! The true and false validators of a Bool Validator Dependency "
    "must be the same type! " <<std::endl << std::endl);
  
}

}

