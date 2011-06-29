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

