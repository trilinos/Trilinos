// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_DependencySheet.hpp"


namespace Teuchos{


DependencySheet::DependencySheet():
  name_("DEP_ANONYMOUS")
{}

DependencySheet::DependencySheet(const std::string &name):
  name_(name)
{}

void DependencySheet::addDependency(RCP<Dependency> dependency){
  Dependency::ConstParameterEntryList dependees = dependency->getDependees();
  for(
    Dependency::ConstParameterEntryList::iterator it  = dependees.begin();
    it != dependees.end();
    ++it)
  {
    dependenciesMap_[*it].insert(dependency);
  }
  dependencies_.insert(dependency);
}

void DependencySheet::removeDependency(RCP<Dependency> dependency){
  Dependency::ConstParameterEntryList dependees = dependency->getDependees();
  for(
    Dependency::ConstParameterEntryList::iterator it  = dependees.begin();
    it != dependees.end();
    ++it)
  {
    dependenciesMap_[*it].erase(dependency);
  }
  dependencies_.erase(dependency);
}

RCP<const DependencySheet::DepSet> DependencySheet::getDependenciesForParameter(
  RCP<const ParameterEntry> dependee) const
{
  if(dependenciesMap_.find(dependee) != dependenciesMap_.end()){
    return rcpFromRef(dependenciesMap_.find(dependee)->second);
  }
  return null;
}

void DependencySheet::printDeps(std::ostream& out) const {
  out << "Dependency Sheet: " << name_ << std::endl << std::endl;
  for(DepSet::const_iterator it = depBegin(); it != depEnd(); ++it){
    (*it)->print(out);
  }
}

void DependencySheet::addDependencies(RCP<DependencySheet> otherSheet){
  for(
    DepSet::iterator it = otherSheet->depBegin();
    it != otherSheet->depEnd();
    ++it
  )
  {
    addDependency(*it);
  }
}


} //namespace Teuchos

