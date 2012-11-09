// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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

