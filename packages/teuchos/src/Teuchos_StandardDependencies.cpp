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
#include "Teuchos_ArrayHelperFunctions.hpp"

namespace Teuchos{

VisualDependency::VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf)
:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::VisualDep), showIf(showIf){}

VisualDependency::VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents, bool showIf)
:Dependency(dependeeName, dependeeParentList, dependents,  Dependency::VisualDep), showIf(showIf){}

VisualDependency::VisualDependency(ParameterParentMap dependees, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf)
:Dependency(dependees, dependentName, dependentParentList, Dependency::VisualDep), showIf(showIf){}

VisualDependency::VisualDependency(ParameterParentMap dependees, ParameterParentMap dependents, bool showIf) 
:Dependency(dependees, dependents, Dependency::VisualDep), showIf(showIf){}

VisualDependency::VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents)
:Dependency(dependeeName, dependeeParentList, dependents,  Dependency::VisualDep){}

VisualDependency::VisualDependency(ParameterParentMap dependees, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList)
:Dependency(dependees, dependentName, dependentParentList, Dependency::VisualDep){}

VisualDependency::VisualDependency(ParameterParentMap dependees, ParameterParentMap dependents) 
:Dependency(dependees, dependents, Dependency::VisualDep){}

bool VisualDependency::isDependentVisible(){
	return dependentVisible;
}

ValidatorDependency::ValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList)
:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::ValidatorDep){}

ValidatorDependency::ValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents)
:Dependency(dependeeName, dependeeParentList, dependents, Dependency::ValidatorDep){}


StringVisualDependency::StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, std::string value, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList, showIf), values(ValueList(1,value)){
	validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, std::string dependentName, 
Teuchos::RCP<Teuchos::ParameterList> parentList, std::string value, bool showIf)
:VisualDependency(dependeeName, parentList, dependentName, parentList, showIf), values(ValueList(1,value)){
	validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, const ValueList& values, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList, showIf), values(values){
	validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, std::string dependentName, 
Teuchos::RCP<Teuchos::ParameterList> parentList, const ValueList& values, bool showIf)
:VisualDependency(dependeeName, parentList, dependentName, parentList, showIf), values(values){
	validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList, 
ParameterParentMap dependents, std::string value, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependents, showIf), values(ValueList(1,value)){
	validateDep();
}

StringVisualDependency::StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents, const ValueList& values, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependents, showIf), values(values){
	validateDep();
}

void StringVisualDependency::evaluate(){
	std::string dependeeValue = getFirstDependeeValue<std::string>();
	ValueList::const_iterator result;
	result = find(values.begin(), values.end(), dependeeValue);
	if((result != values.end() && showIf) || (result==values.end() && !showIf)){
		dependentVisible = true;
	}
	else{
		dependentVisible = false;
	}
}

void StringVisualDependency::validateDep(){
	/*
	 * This error should never get thrown, unless someone
	 * is doing something wonky in a sublcass.
	 */
	if(dependees.size() != 1){
		throw InvalidDependencyException("Uh oh. Looks like you tried to make a "
		"String Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
		"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
		"be doing some subclassing you sly-dog ;)\n\n" 
		"Error: A String Visual Dependency must have exactly 1 dependee. " 
		"You have tried to assign it "+ QString::number(dependees.size()).toStdString() + " dependees.\n" 
		"Dependees: " + getDependeeNamesString() + "\n" 
		"Dependents: " + getDependentNamesString());
	}
	if(!getFirstDependee()->isType<std::string>()){
		throw InvalidDependencyException("Ay no! The dependee of a "
		"String Visual Dependency must be of type string.\n"
		"Problem dependee: " + getFirstDependeeName() + "\n"
		"Actual type: " + getFirstDependee()->getAny().typeName() + "\n"
		"Dependent: " + getDependentNamesString());
	}
}

BoolVisualDependency::BoolVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList, showIf){
	validateDep();
}

BoolVisualDependency::BoolVisualDependency(std::string dependeeName, std::string dependentName, 
Teuchos::RCP<Teuchos::ParameterList> parentList, bool showIf)
:VisualDependency(dependeeName, parentList, dependentName, parentList, showIf){
	validateDep();
}

BoolVisualDependency::BoolVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents, bool showIf)
:VisualDependency(dependeeName, dependeeParentList, dependents, showIf){
	validateDep();
}

void BoolVisualDependency::evaluate(){
	bool dependeeValue = getFirstDependeeValue<bool>();
	if((dependeeValue && showIf) || (!dependeeValue && !showIf)){
		dependentVisible = true;
	}
	else{
		dependentVisible = false;
	}
}

void BoolVisualDependency::validateDep(){
	/*
	 * This error should never get thrown, unless someone
	 * is doing something wonky in a sublcass.
	 */
	if(dependees.size() != 1){
		throw InvalidDependencyException("Uh oh. Looks like you tried to make a "
		"Bool Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
		"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
		"be doing some subclassing you sly-dog ;)\n\n" 
		"Error: A Bool Visual Dependency must have exactly 1 dependee. " 
		"You have tried to assign it "+ QString::number(dependees.size()).toStdString() + " dependees.\n" 
		"Dependees: " + getDependeeNamesString() + "\n" 
		"Dependents: " + getDependentNamesString());
	}

	if(!getFirstDependee()->isType<bool>()){
		throw InvalidDependencyException("Ay no! The dependee of a "
		"Bool Visual Dependency must be of type bool.\n"
		"Problem dependee: " + getFirstDependeeName() + "\n"
		"Actual type: " + getFirstDependee()->getAny().typeName() + "\n"
		"Dependents: " + getDependentNamesString());
	}
}

ConditionVisualDependency::ConditionVisualDependency(Teuchos::RCP<Condition> condition,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf)
:VisualDependency(condition->getAllParameters(), dependentName, dependentParentList, showIf), condition(condition){
	validateDep();
}

ConditionVisualDependency::ConditionVisualDependency(Teuchos::RCP<Condition> condition, ParameterParentMap dependents, bool showIf)
:VisualDependency(condition->getAllParameters(), dependents, showIf), condition(condition){
	validateDep();
}

void ConditionVisualDependency::evaluate(){
	bool conditionState = condition->isConditionTrue();
	if((conditionState && showIf) || (!conditionState && !showIf)){
		dependentVisible = true;
	}
	else{
		dependentVisible = false;
	}
}

void ConditionVisualDependency::validateDep(){} 

NumberArrayLengthDependency::NumberArrayLengthDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, int (*func)(int))
:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberArrayLengthDep), func(func){
	validateDep();
}

NumberArrayLengthDependency::NumberArrayLengthDependency(std::string dependeeName, std::string dependentName, 
Teuchos::RCP<Teuchos::ParameterList> parentList, int (*func)(int))
:Dependency(dependeeName, parentList, dependentName, parentList, Dependency::NumberArrayLengthDep), func(func){
	validateDep();
}

NumberArrayLengthDependency::NumberArrayLengthDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents, int (*func)(int))
:Dependency(dependeeName, dependeeParentList, dependents, Dependency::NumberArrayLengthDep), func(func){
	validateDep();
}

int NumberArrayLengthDependency::runFunction(int argument) const{
	if(func !=0){
		return (*func)(argument);
	}
	else{
		return argument;
	}
}

template <class S>
void NumberArrayLengthDependency::modifyArrayLength(int newLength, Teuchos::ParameterEntry* dependentToModify){
	const Teuchos::Array<S> originalArray = Teuchos::any_cast<Teuchos::Array<S> >(dependentToModify->getAny()); 
	Teuchos::RCP<const Teuchos::EnhancedNumberValidator<S> > potentialValidator = Teuchos::null;
	if(!Teuchos::is_null(dependentToModify->validator())){
			potentialValidator = Teuchos::rcp_dynamic_cast<const Teuchos::ArrayNumberValidator<S> >(dependentToModify->validator(),true)->getPrototype();
	}
	Teuchos::Array<S> newArray;
	int i;
	for(i=0; i<originalArray.size() && i<newLength; ++i){
		newArray.append(originalArray[i]);
	}
	for(;i<newLength;++i){
		if(Teuchos::is_null(potentialValidator)){
			newArray.append(0);
		}
		else{
			newArray.append(potentialValidator->getMin());
		}
	}
	dependentToModify->setValue(newArray, false, dependentToModify->docString(), dependentToModify->validator());
}

template<>
void NumberArrayLengthDependency::modifyArrayLength<std::string>(int newLength, Teuchos::ParameterEntry* dependentToModify){
	const Teuchos::Array<std::string> originalArray = Teuchos::any_cast<Teuchos::Array<std::string> >(dependentToModify->getAny()); 
	Teuchos::Array<std::string> newArray;
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = dependentToModify->validator();
	int i;
	for(i=0; i<originalArray.size() && i<newLength; ++i){
		newArray.append(originalArray[i]);
	}
	for(;i<newLength;++i){
		if(Teuchos::is_null(validator)){
			newArray.append(" ");
		}
		else if(!Teuchos::is_null(Teuchos::rcp_dynamic_cast<const Teuchos::ArrayFileNameValidator>(validator))){
			newArray.append(QDir::homePath().toStdString());
		}
		else{
			newArray.append(validator->validStringValues()->at(0));
		}
	}
	dependentToModify->setValue(newArray, false, dependentToModify->docString(), dependentToModify->validator());
}

void NumberArrayLengthDependency::evaluate(){
	int newLength = runFunction(getFirstDependeeValue<int>());
	if(newLength <0){
		std::stringstream oss;
		std::string msg;
		oss << "Ruh Roh Shaggy! Looks like a dependency tried to set the length of the Array(s) in the " <<
		getDependentNamesString() << " parameter(s) to a negative number. Silly. You can't have an Array " <<
		"with a negative length! You should probably contact the maintainer of this program and " <<
		"give him or her the following information: \n\n" <<
		"Error:\n" <<
		"An attempt was made to set the length of an Array to a negative number by a NumberArrayLengthDependency\n" <<
		"Dependency Type: " << QString::number(getType()).toStdString() + "\n" <<
		"Problem Dependee: " << getFirstDependeeName() <<
		"Problem Dependents: " << getDependentNamesString();
		msg = oss.str();
		throw Teuchos::Exceptions::InvalidParameterValue(msg);
	}
	Teuchos::ParameterEntry *currentDependent;
	for(ParameterParentMap::const_iterator it = dependents.begin(); it != dependents.end(); ++it){ 
		currentDependent = it->second->getEntryPtr(it->first);
		//QString dependentType = determineArrayType(currentDependent);
		if(currentDependent->getAny().type() == typeid(Teuchos::Array<int>)){
			modifyArrayLength<int>(newLength, currentDependent);
		}
		else if(currentDependent->getAny().type() == typeid(Teuchos::Array<short>)){
			modifyArrayLength<short>(newLength, currentDependent);
		}
		else if(currentDependent->getAny().type() == typeid(Teuchos::Array<double>)){
			modifyArrayLength<double>(newLength, currentDependent);
		}
		else if(currentDependent->getAny().type() == typeid(Teuchos::Array<float>)){
			modifyArrayLength<float>(newLength, currentDependent);
		}
		else if(currentDependent->getAny().type() == typeid(Teuchos::Array<std::string>)){
			modifyArrayLength<std::string>(newLength, currentDependent);
		}
	}
}

void NumberArrayLengthDependency::validateDep(){
	/*
	 * This error should never get thrown, unless someone
	 * is doing something wonky in a sublcass.
	 */
	if(dependees.size() != 1){
		throw InvalidDependencyException("Uh oh. Looks like you tried to make a "
		"Number Array Length Dependency doesn't have exactly one currentDependee. This is kind of a problem. " 
		"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
		"be doing some subclassing you sly-dog ;)\n\n" 
		"Error: A Number Array Length Dependency must have exactly 1 currentDependee. " 
		"You have tried to assign it "+ QString::number(dependees.size()).toStdString() + " dependees.\n" 
		"Dependees: " + getDependeeNamesString() + "\n" 
		"Dependents: " + getDependentNamesString());
	}

	if(!getFirstDependee()->isType<int>() && !getFirstDependee()->isType<short>()){
		throw InvalidDependencyException("Ay no! The currentDependee for an "
		"Array Length Dependency must be either of type short or of type int.\n" 
		"Problem currentDependee: " + getFirstDependeeName() + "\n" 
		"Actual type: " + getFirstDependee()->getAny().typeName() + "\n"
		"Dependent: " + getDependentNamesString());
	}
	Teuchos::ParameterEntry *currentDependent;
	for(ParameterParentMap::const_iterator it = dependents.begin(); it != dependents.end(); ++it){ 
		currentDependent = it->second->getEntryPtr(it->first);
		if(!doesParameterContainArray(currentDependent)){
			throw InvalidDependencyException("Ay no! The dependent of an "
			"Array Length Dependency must be an array.\n"
			"Problem dependent: " + it->first + "\n" 
			"Actual type: " + currentDependent->getAny().typeName() + "\n"
			"Dependees: " + getDependeeNamesString());
		}
	}
}

StringValidatorDependency::StringValidatorDependency(std::string currentDependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,
ValueToValidatorMap valuesAndValidators, Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
:ValidatorDependency(currentDependeeName, dependeeParentList, dependentName, dependentParentList), defaultValidator(defaultValidator), valuesAndValidators(valuesAndValidators){
	validateDep();
}

StringValidatorDependency::StringValidatorDependency(std::string currentDependeeName, std::string dependentName, 
Teuchos::RCP<Teuchos::ParameterList> parentList, ValueToValidatorMap valuesAndValidators, 
Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
:ValidatorDependency(currentDependeeName, parentList, dependentName, parentList), defaultValidator(defaultValidator), valuesAndValidators(valuesAndValidators){
	validateDep();
}

StringValidatorDependency::StringValidatorDependency(std::string currentDependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents,
ValueToValidatorMap valuesAndValidators, Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
:ValidatorDependency(currentDependeeName, dependeeParentList, dependents), defaultValidator(defaultValidator), valuesAndValidators(valuesAndValidators){
	validateDep();
}

void StringValidatorDependency::evaluate(){
	std::string currentDependeeValue = getFirstDependeeValue<std::string>();
	Teuchos::ParameterEntry *currentDependent;
	for(ParameterParentMap::const_iterator it = dependents.begin(); it != dependents.end(); ++it){ 
		currentDependent = it->second->getEntryPtr(it->first);
		if(valuesAndValidators.find(currentDependeeValue) == valuesAndValidators.end()){
			currentDependent->setValidator(defaultValidator);
		}
		else{
			currentDependent->setValidator(valuesAndValidators[currentDependeeValue]);
		}
	}
}
void StringValidatorDependency::validateDep(){
	/*
	 * This error should never get thrown, unless someone
	 * is doing something wonky in a sublcass.
	 */
	if(dependees.size() != 1){
		throw InvalidDependencyException("Uh oh. Looks like you tried to make a "
		"String Validator Dependency doesn't have exactly one dependee. This is kind of a problem. " 
		"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
		"be doing some subclassing you sly-dog ;)\n\n" 
		"Error: A String Validator Dependency must have exactly 1 dependee. "
		"You have tried to assign it "+ QString::number(dependees.size()).toStdString() + " dependees.\n"
		"Dependees: " + getDependeeNamesString() + "\n"
		"Dependents: " + getDependentNamesString());
	}

	if(!getFirstDependee()->isType<std::string>()){
		throw InvalidDependencyException("Ay no! The dependee of a "
		"String Validator Dependency must be of type string.\n"
		"Problem dependee: " + getFirstDependeeName() + "\n"
		"Actual type: " + getFirstDependee()->getAny().typeName() + "\n"
		"Dependent: " + getDependentNamesString());
	}
	for(ValueToValidatorMap::const_iterator it = valuesAndValidators.begin(); it != valuesAndValidators.end(); ++it){
		Teuchos::ParameterEntry *currentDependent;
		for(ParameterParentMap::const_iterator it2 = dependents.begin(); it2 != dependents.end(); ++it2){ 
			currentDependent = it2->second->getEntryPtr(it2->first);
			if(typeid(*(currentDependent->validator().get())) != typeid(*(it->second.get()))){
				throw InvalidDependencyException("Ay no! The validator of a dependent of a "
				"String Validator Dependency must be the same type as all of the validators "
				"in the valuesAndValidators map. "
				"Note this means that the dependent must have an initial validator.\n"
				"Problem dependent: " + it->first + "\n"
				"Validator Type: " + typeid(*(currentDependent->validator())).name() + "\n"
				"Problem Map Value: " + it->first + "\n"
				"The Validator type associated with the Map Value: " + typeid(*(it->second)).name());
			}
		}
	}
}

BoolValidatorDependency::BoolValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,
Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator, Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator)
:ValidatorDependency(dependeeName, dependeeParentList, dependentName, dependentParentList), trueValidator(trueValidator), falseValidator(falseValidator){
	validateDep();
}

BoolValidatorDependency::BoolValidatorDependency(std::string dependeeName, std::string dependentName, 
Teuchos::RCP<Teuchos::ParameterList> parentList, Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator,
Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator)
:ValidatorDependency(dependeeName, parentList, dependentName, parentList), trueValidator(trueValidator), falseValidator(falseValidator){
	validateDep();
}

BoolValidatorDependency::BoolValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
ParameterParentMap dependents,
Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator, Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator)
:ValidatorDependency(dependeeName, dependeeParentList, dependents), trueValidator(trueValidator), falseValidator(falseValidator){
	validateDep();
}

void BoolValidatorDependency::evaluate(){
	bool dependeeValue = getFirstDependeeValue<bool>();
	Teuchos::ParameterEntry *currentDependent;
	for(ParameterParentMap::const_iterator it = dependents.begin(); it != dependents.end(); ++it){ 
		currentDependent = it->second->getEntryPtr(it->first);
		dependeeValue ? currentDependent->setValidator(trueValidator) : currentDependent->setValidator(falseValidator);
	}
}

void BoolValidatorDependency::validateDep(){
	/*
	 * This error should never get thrown, unless someone
	 * is doing something wonky in a sublcass.
	 */
	if(dependees.size() != 1){
		throw InvalidDependencyException("Uh oh. Looks like you tried to make a "
		"Bool Validator Dependency doesn't have exactly one dependee. This is kind of a problem. " 
		"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
		"be doing some subclassing you sly-dog ;)\n\n" 
		"Error: A Bool Validator Dependency must have exactly 1 dependee. "
		"You have tried to assign it "+ QString::number(dependees.size()).toStdString() + " dependees.\n"
		"Dependees: " + getDependeeNamesString() + "\n" 
		"Dependents: " + getDependentNamesString());
	}

	if(!getFirstDependee()->isType<bool>()){
		throw InvalidDependencyException("Ay no! The dependee of a "
		"Bool Validator Dependency must be of type boolean.\n"
		"Problem dependee: " + getFirstDependeeName() + "\n"
		"Actual type: " + getFirstDependee()->getAny().typeName() + "\n"
		"Dependent: " + getDependentNamesString());
	}
	Teuchos::ParameterEntry *currentDependent;
	for(ParameterParentMap::const_iterator it = dependents.begin(); it != dependents.end(); ++it){ 
		currentDependent = it->second->getEntryPtr(it->first);
		if(typeid(*(currentDependent->validator().get())) != typeid(*(trueValidator.get()))){
			throw InvalidDependencyException("Ay no! The validator of a dependent of a "
			"Bool Validator Dependency must be the same type as the \"true\" validator. "
			"Note this means that the dependent must have an initial validator.\n"
			"Problem dependent: " + it->first + "\n"
			"Validator Type: " + typeid(*(currentDependent->validator().get())).name() + "\n"
			"Type of the \"true\" validator: " + typeid(*(trueValidator.get())).name());
		}
		if(typeid(*(currentDependent->validator().get())) != typeid(*(falseValidator.get()))){
			throw InvalidDependencyException("Ay no! The validator of a dependent of a "
			"Bool Validator Dependency must be the same type as the \"false\" validator. "
			"Note this means that the dependent must have an initial validator.\n"
			"Problem dependent: " + it->first + "\n"
			"Validator Type: " + typeid(*(currentDependent->validator().get())).name() + "\n"
			"Type of the \"false\" validator: " + typeid(*(falseValidator.get())).name());
		}
	}
}

}

