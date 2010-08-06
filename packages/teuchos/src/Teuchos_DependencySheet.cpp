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


DependencySheet::DependencySheet(RCP<ParameterList> rootList):
  name_("DEP_ANONYMOUS")

DependencySheet::DependencySheet(const std::string &name):
  name_(name)

bool DependencySheet::addDependency(RCP<Dependency> dependency){
  bool successfulAdd = true;
  Dependency::ParameterParentMap dependees = dependency->getDependees();
  for(Dependency::ParameterParentMap::iterator it  = dependees.begin(); it != dependees.end(); ++it){
    if(!dependencies_[it->second->getEntryPtr(it->first)].insert(dependency).second){
      successfulAdd = false;
    }
  }
  return successfulAdd;
}

bool DependencySheet::removeDependency(RCP<Dependency> dependency){
  bool succesfulRemove = true;  
  Dependency::ParameterParentMap dependees = dependency->getDependees();
  for(Dependency::ParameterParentMap::iterator it  = dependees.begin(); it != dependees.end(); ++it){
    bool successFullCurrentRemove = false;
    ParameterEntry* currentDependee = it->second->getEntryPtr(it->first);
    for(
      DepSet::iterator it2 = dependencies_[currentDependee].begin(); 
      it2 != dependencies_[currentDependee].end(); 
    ++it2)
    {
      if((*it2) == dependency){
        dependencies_[currentDependee].erase(it2);
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

void DependencySheet::printDeps(){
  std::cout << "Dependency Sheet: " << name_ << "\n\n";
  for(DepMap::iterator it = depBegin(); it != depEnd(); ++it){
    const ParameterEntry* dependee = it->first;
    for(DepSet::iterator it2 = dependencies_.find(dependee)->second.begin(); it2 != dependencies_.find(dependee)->second.end(); ++it2){
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


} //namespace Teuchos

