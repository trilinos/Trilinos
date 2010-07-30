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



#include "Teuchos_Dependency.hpp"

namespace Teuchos{

Dependency::Dependency(ParameterParentMap& dependees, ParameterParentMap& dependents, Type type):
	type(type)
{
	intitializeDependeesAndDependents(dependees,dependents);
}

Dependency::Dependency(ParameterParentMap& dependees, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, Type type):
	type(type)
{
	ParameterParentMap dependents;
	dependents.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(dependentName, dependentParentList));
	intitializeDependeesAndDependents(dependees,dependents);
}

Dependency::Dependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap& dependents, Type type):
	type(type)
{
	ParameterParentMap dependees;
	dependees.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(dependeeName, dependeeParentList));
	intitializeDependeesAndDependents(dependees,dependents);
}
	
Dependency::Dependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
 	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, Type type):
	type(type)
{
	ParameterParentMap dependees;
	dependees.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(dependeeName, dependeeParentList));
	ParameterParentMap dependents;
	dependents.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(dependentName, dependentParentList));
	intitializeDependeesAndDependents(dependees, dependents);
}

void Dependency::intitializeDependeesAndDependents(ParameterParentMap& dependees, ParameterParentMap& dependents){
	ParameterParentMap::iterator it;
	for(it = dependees.begin(); it != dependees.end(); ++it){
		if(it->second->getEntryPtr(it->first) == NULL){
			throw InvalidDependencyException("The Dependee Parameter \"" + it->first + "\" does "
			"not exist in the given Dependent Parent List \"" + it->second->name() + "\"."
			"\n\nBummer! Maybe you just mispelled something? Why not go back and check to make sure "
			"you've got the names of the dependee and the depedent right? "
			"It might also be that you just didn't specify the correct parent lists for the dependent and "
			"dependee. Either way, I'm sure it's just a simple mistake. "
			"You're a great programmer, I'm sure you'll figure it out! :)");
		}
		else{
			this->dependees.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(it->first, it->second));
			dependeeNames.insert(it->first);
		}
	}
	for(it = dependents.begin(); it != dependents.end(); ++it){
		if(it->second->getEntryPtr(it->first) == NULL){
			throw InvalidDependencyException("The Dependent Parameter \"" + it->first + "\" does "
			"not exist in the given Dependent Parent List \"" + it->second->name() + "\"."
			"\n\nBummer! Maybe you just mispelled something? Why not go back and check to make sure "
			"you've got the names of the dependee and the depedent right? "
			"It might also be that you just didn't specify the correct parent lists for the dependent and "
			"dependee. Either way, I'm sure it's just a simple mistake. "
			"You're a great programmer, I'm sure you'll figure it out! :)");
		}
		else{
			this->dependents.insert(std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> >(it->first, it->second));
			dependentNames.insert(it->first);
		}
	}
}


Dependency::ParameterParentMap Dependency::getDependees() const{
	return dependees;
}

Dependency::ParameterParentMap Dependency::getDependents() const{
	return dependents;
}

std::set<std::string> Dependency::getDependeeNames() const{
	return dependeeNames;
}

std::string Dependency::getDependeeName(const Teuchos::ParameterEntry* dependee) const{
	for(ParameterParentMap::const_iterator it = dependees.begin(); it != dependees.end(); ++it){
		if(dependee == it->second->getEntryPtr(it->first)){
			return it->first;
		}
	}
	throw InvalidDependencyException("Fooey! Looks like you tried to get the name of a dependee parameter "
	"that isn't actually a dependee parameter for this dependency. Make sure you're giving this funciton "
	"the right pointer. Maybe the information below can help you out.\n\n"
	"Error: Dependency does not contain specified dependee parameter.\n"
	"Dependee(s): " + getDependeeNamesString() + "\n"
	"Dependent(s): " + getDependentNamesString() + "\n");
	return "";
}


std::set<std::string> Dependency::getDependentNames() const{
	return dependentNames;
}

std::string Dependency::getDependeeNamesString() const{
	std::string names = "";
	for(std::set<std::string>::const_iterator it=dependeeNames.begin(); it != dependeeNames.end(); ++it){
		names += *it + " ";
	}
	return names;
}

std::string Dependency::getDependentNamesString() const{
	std::string names = "";
	for(std::set<std::string>::const_iterator it=dependentNames.begin(); it != dependentNames.end(); ++it){
		names += *it + " ";
	}
	return names;
}

Dependency::Type Dependency::getType() const{
	return type;
}

bool Dependency::doesListContainList(Teuchos::RCP<Teuchos::ParameterList> parentList, Teuchos::RCP<Teuchos::ParameterList> listToFind){
	if(parentList.get() == listToFind.get()){
		return true;
	}
	else{
		for(Teuchos::ParameterList::ConstIterator it = parentList->begin(); it!=parentList->end(); ++it){
			if(it->second.isList()){
				if(doesListContainList(sublist(parentList, it->first,true), listToFind)){
					return true;
				}
			}
		}
	}
	return false;
}

}

