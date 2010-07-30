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



#include "Teuchos_StandardConditions.hpp"

namespace Teuchos{

ParameterCondition::ParameterCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, bool whenParamEqualsValue):
	Condition(Condition::ParamCon),
	parameterName(parameterName),
	parentList(parentList),
	whenParamEqualsValue(whenParamEqualsValue)
{
	parameter = parentList->getEntryPtr(parameterName);
	if(parameter == NULL){
		throw InvalidConditionException("Oh noes!!!!!!! Looks like the parameter " +
		parameterName + " isn't actually contained in the " + parentList->name() + " parameter list. "
		"You should go back and check your code. Maybe the information below can help you.\n\n"
		"Error: Parameter not contained in specified Parent list:\n"
		"Problem Parameter: " + parameterName + "\n"
		"Problem List: " + parentList->name());
	}
}

bool ParameterCondition::containsAtLeasteOneParameter(){
	return true;
}

Dependency::ParameterParentMap ParameterCondition::getAllParameters(){
	Dependency::ParameterParentMap toReturn;
	toReturn.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(parameterName, parentList));
	return toReturn;
}

BinaryLogicalCondition::BinaryLogicalCondition(ConditionList& conditions):
	Condition(Condition::BinLogicCon),
	conditions(conditions)
{
	if(conditions.size() ==0){
		throw InvalidConditionException("Sorry bud, but you gotta at least give me one condition "
		"when you're constructing a BinaryLogicalCondition. Looks like you didn't. I'm just gonna "
		"chalk it up a silly little mistake though. Take a look over your conditions again and make sure "
		"you don't ever give any of your BinaryLogicConditions and empty condition list.\n\n"
		"Error: Empty condition list given to a BinaryLogicalCondition constructor.");
	}
}


void BinaryLogicalCondition::addCondition(Teuchos::RCP<Condition> toAdd){
	conditions.append(toAdd);
}

bool BinaryLogicalCondition::containsAtLeasteOneParameter(){
	for(ConditionList::iterator it=conditions.begin(); it!=conditions.end(); ++it){
		if((*it)->containsAtLeasteOneParameter()){
			return true;
		}
	}
	return false;
}

Dependency::ParameterParentMap BinaryLogicalCondition::getAllParameters(){
	Dependency::ParameterParentMap toReturn;
	Dependency::ParameterParentMap currentMap;
	for(ConditionList::iterator it = conditions.begin(); it != conditions.end(); ++it){
		currentMap = (*it)->getAllParameters();
		toReturn.insert(currentMap.begin(), currentMap.end());
	}
	return toReturn;
}

OrCondition::OrCondition(ConditionList& conditions):
	BinaryLogicalCondition(conditions){}

bool OrCondition::isConditionTrue(){
	ConditionList::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn || (*it)->isConditionTrue());
	}
	return toReturn;
}

AndCondition::AndCondition(ConditionList& conditions):
	BinaryLogicalCondition(conditions){}

bool AndCondition::isConditionTrue(){
	ConditionList::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn && (*it)->isConditionTrue());
	}
	return toReturn;
}

EqualsCondition::EqualsCondition(ConditionList& conditions):
	BinaryLogicalCondition(conditions){}

bool EqualsCondition::isConditionTrue(){
	ConditionList::iterator it = conditions.begin();
	bool toReturn = (*it)->isConditionTrue();
	++it;
	for(;it != conditions.end(); ++it){
		toReturn = (toReturn == (*it)->isConditionTrue());
	}
	return toReturn;
}

NotCondition::NotCondition(Teuchos::RCP<Condition> condition):
	Condition(Condition::NotCon),
	condition(condition)
{
	if(condition.is_null()){
		throw InvalidConditionException("OOOOOOOOPppppps! Looks like you tried to give me "
		"a null pointer when you were making a not conditon. That's a no no. Go back and "
		"checkout your not conditions and make sure you didn't give any of them a null pointer "
		"as an argument to the constructor.\n\n"
		"Error: Null pointer given to NotCondition constructor.");
	}
}

bool NotCondition::isConditionTrue(){
	return (!condition->isConditionTrue());
}

bool NotCondition::containsAtLeasteOneParameter(){
	return condition->getType() == Condition::ParamCon;
}

Dependency::ParameterParentMap NotCondition::getAllParameters(){
	return condition->getAllParameters();
}

StringCondition::StringCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, std::string value, bool whenParamEqualsValue):
	ParameterCondition(parameterName, parentList, whenParamEqualsValue), values(ValueList(1,value))
{
	if(!parameter->isType<std::string>()){
		throw InvalidConditionException("The parameter of a String Condition "
		"must be of type string. \n"
		"Expected type: std::string\n"
		"Actual type: " + parameter->getAny().typeName() + "\n");
	}
}

StringCondition::StringCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, ValueList values, bool whenParamEqualsValue):
	ParameterCondition(parameterName, parentList, whenParamEqualsValue), values(values)
{
	if(!parameter->isType<std::string>()){
		throw InvalidConditionException("The parameter of a String Condition "
		"must be of type string. \n"
		"Expected type: std::string\n"
		"Actual type: " + parameter->getAny().typeName() + "\n");
	}
}

bool StringCondition::isConditionTrue(){
	bool toReturn = find(values.begin(), values.end(), Teuchos::getValue<std::string>(*parameter)) != values.end();
	return whenParamEqualsValue ? toReturn : !toReturn;
}

BoolCondition::BoolCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, bool whenParamEqualsValue):
	ParameterCondition(parameterName, parentList, whenParamEqualsValue)
{
	if(!parameter->isType<bool>()){
		throw InvalidConditionException("The parameter of a Bool Condition "
		"must be of type Bool. \n"
		"Expected type: Bool\n"
		"Actual type: " + parameter->getAny().typeName() + "\n");
	}
}

bool BoolCondition::isConditionTrue(){
	bool toReturn = Teuchos::getValue<bool>(*parameter);
	return whenParamEqualsValue ? toReturn : !toReturn;
}

}

