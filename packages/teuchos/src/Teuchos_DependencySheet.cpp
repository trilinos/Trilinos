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



#include "Teuchos_DependencySheet.hpp"

namespace Teuchos{


DependencySheet::DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList):name("DEP_ANONYMOUS"), rootList(rootList){}

DependencySheet::DependencySheet(Teuchos::RCP<Teuchos::ParameterList> rootList, const std::string &name):name(name), rootList(rootList){}

bool DependencySheet::addDependency(Teuchos::RCP<Teuchos::Dependency> dependency){
	validateExistanceInRoot(dependency);
	bool successfulAdd = true;
	Dependency::ParameterParentMap dependees = dependency->getDependees();
	for(Dependency::ParameterParentMap::iterator it  = dependees.begin(); it != dependees.end(); ++it){
		if(!dependencies[it->second->getEntryPtr(it->first)].insert(dependency).second){
			successfulAdd = false;
		}
	}
	return successfulAdd;
}

bool DependencySheet::removeDependency(Teuchos::RCP<Teuchos::Dependency> dependency){
	bool succesfulRemove = true;	
	Dependency::ParameterParentMap dependees = dependency->getDependees();
	for(Dependency::ParameterParentMap::iterator it  = dependees.begin(); it != dependees.end(); ++it){
		bool successFullCurrentRemove = false;
		Teuchos::ParameterEntry* currentDependee = it->second->getEntryPtr(it->first);
		for(DepSet::iterator it2 = dependencies[currentDependee].begin(); it2 != dependencies[currentDependee].end(); ++it2){
			if((*it2) == dependency){
				dependencies[currentDependee].erase(it2);
				successFullCurrentRemove = true;
				break;
			}
		}
		if(!successFullCurrentRemove){
			succesfulRemove = false;
		}
	}
	return succesfulRemove;
}

bool DependencySheet::hasDependents(const Teuchos::ParameterEntry* dependee) const{
	return (dependencies.find(dependee) != dependencies.end() && dependencies.find(dependee)->second.size() > 0);
}

const DependencySheet::DepSet& DependencySheet::getDependenciesForParameter(const Teuchos::ParameterEntry* dependee) const{
	return dependencies.find(dependee)->second;
}

DependencySheet::DepMap::iterator DependencySheet::depBegin(){
	return dependencies.begin();
}

DependencySheet::DepMap::iterator DependencySheet::depEnd(){
	return dependencies.end();
}

DependencySheet::DepMap::const_iterator DependencySheet::depBegin() const{
	return dependencies.begin();
}

DependencySheet::DepMap::const_iterator DependencySheet::depEnd() const{
	return dependencies.end();
}

void DependencySheet::printDeps(){
	std::cout << "Dependency Sheet: " << name << "\n\n";
	for(DepMap::iterator it = depBegin(); it != depEnd(); ++it){
		const Teuchos::ParameterEntry* dependee = it->first;
		for(DepSet::iterator it2 = dependencies.find(dependee)->second.begin(); it2 != dependencies.find(dependee)->second.end(); ++it2){
			std::set<std::string> dependeeNames = (*it2)->getDependeeNames();
			std::cout << "Dependees: \n";
			for(std::set<std::string>::iterator it3 = dependeeNames.begin(); it3 != dependeeNames.end(); ++it3){
				std::cout << "\t" << *it3 << "\n";
			}
			std::set<std::string> dependentNames = (*it2)->getDependentNames();
			std::cout << "Dependents: \n";
			for(std::set<std::string>::iterator it3 = dependentNames.begin(); it3 != dependentNames.end(); ++it3){
				std::cout << "\t" << *it3 << "\n";
			}
			std::cout << "\n";
			std::cout << "Type: " << (*it2)->getType() << "\n\n";
		}
	}
}

void DependencySheet::validateExistanceInRoot(Teuchos::RCP<Teuchos::Dependency> dependency){
	Dependency::ParameterParentMap::const_iterator it;
	Teuchos::ParameterEntry *currentDependee;
	Dependency::ParameterParentMap dependees = dependency->getDependees();
	for(it = dependees.begin(); it != dependees.end(); ++it){ 
		currentDependee = it->second->getEntryPtr(it->first);
		TEST_FOR_EXCEPTION(!Dependency::doesListContainList(rootList, it->second),
			InvalidDependencyException,
			"FAILED TO ADD DEPENDENCY!\n\n"
			"Sorry for the yelling there, but this is kind of a big deal. Dependencies are hard and complex so don't beat "
			"yourself up too much. Mistakes are easy to make when dealing with dependencies. "
			"And besides, I'm gonna do my best to help you out! I'm sure with the informationg below you'll be able to figure out what "
			"exactly went wrong. I've got confidence in you! :)\n\n"
			"Error:\n"
			"An attempt was made to add a dependency containing a the dependee parameter \"" << it->first << "\""
			" to the Dependency Sheet \"" << name << "\"."
			" The Dependency Sheet's root list does not contain nor does it have"
			" child ParameterLists that contain the parameter.\n"
			"Dependency Sheet: " << name << "\n"
			"Dependency Type: " <<dependency->getType() << "\n"
			"Bad Dependee Name: " << it->first);
		}
	}

	Teuchos::ParameterEntry *currentDependent;
	Dependency::ParameterParentMap dependents = dependency->getDependents();
	for(it = dependents.begin(); it != dependents.end(); ++it){ 
		currentDependent = it->second->getEntryPtr(it->first);
		TEST_FOR_EXCEPTION(!Dependency::doesListContainList(rootList, it->second),
			InvalidDependencyException,
			"FAILED TO ADD DEPENDENCY!\n\n"
			"Sorry for the yelling there, but this is kind of a big deal. Dependencies are hard and complex so don't beat "
			"yourself up too much. Mistakes are easy to make when dealing with dependencies. "
			"And besides, I'm gonna do my best to help you out! I'm sure with the informationg below you'll be able to figure out what "
			"exactly went wrong. I've got confidence in you! :)\n\n"
			"Error:\n"
			"An attempt was made to add a dependency containing a the dependent parameter \"" << it->first << "\""
			" to the Dependency Sheet \"" << name << "\"."
			" The Dependency Sheet's root list does not contain nor does it have"
			" child ParameterLists that contain the parameter.\n"
			"Dependency Sheet: " << name << "\n"
			"Dependency Type: " << dependency->getType() << "\n"
			"Bad Dependent Name: " << it->first);
		}
	}
}

}

