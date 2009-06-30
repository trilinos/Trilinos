#include "TivaBuena_Dependency.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TivaBuena_InvalidDependencyException.hpp"

namespace TivaBuena{

Dependency::Dependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, Type type){
	this->dependentParentList = dependentParentList;
	this->dependeeParentList = dependeeParentList;
	this->dependeeName = dependeeName;
	this->dependentName = dependentName;
	if(dependeeParentList->getEntryPtr(dependeeName) == NULL)
		throw InvalidDependencyException("The Dependee Parameter \"" + dependeeName + "\" does "
		"not exist in the given Dependee Parent List \"" + dependeeParentList->name() + "\".");
	else
		this->dependee = dependeeParentList->getEntryPtr(dependeeName);
	if(dependentParentList->getEntryPtr(dependentName) == NULL)
		throw InvalidDependencyException("The Dependent Parameter \"" + dependentName + "\" does "
		"not exist in the given Dependent Parent List \"" + dependentParentList->name() + "\".");
	else
		this->dependent = dependentParentList->getEntryPtr(dependentName);
	this->type = type;
}

bool Dependency::isDependentParentInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList){
	return doesListContainList(potentialParentList, dependentParentList->name());
}

bool Dependency::isDependeeParentInList(Teuchos::RCP<Teuchos::ParameterList> potentialParentList){
	return doesListContainList(potentialParentList, dependeeParentList->name());
}

bool Dependency::doesListContainList(Teuchos::RCP<Teuchos::ParameterList> parentList, std::string listName){
	if(parentList->name() == listName)
		return true;
	else if(parentList->isSublist(listName))
		return true;
	else{
		for(Teuchos::ParameterList::ConstIterator it = parentList->begin(); it!=parentList->end(); it++){
			if(it->second.isList()){
				if(doesListContainList(sublist(parentList, it->first, true), listName))
					return true;
			}
		}
	}
	return false;
}

const Teuchos::ParameterEntry* Dependency::getDependee() const{
	return dependee;
}

const Teuchos::ParameterEntry* Dependency::getDependent() const{
	return dependent;
}


const std::string& Dependency::getDependeeName() const{
	return dependeeName;
}

const std::string& Dependency::getDependentName() const{
	return dependentName;
}
/*
Teuchos::RCP<const Teuchos::ParameterList> Dependency::getDependentParentList() const{
	return dependentParentList;
}

Teuchos::RCP<const Teuchos::ParameterList> Dependency::getDependeeParentList() const{
	return dependeeParentList;
}*/

Dependency::Type Dependency::getType() const{
	return type;
}

}

