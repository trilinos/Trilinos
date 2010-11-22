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


Dependency::Dependency(
ConstParameterEntryList dependees,
ParameterEntryList dependents):
  dependees_(dependees), dependents_(dependents)
{
  checkDependeesAndDependents();
  createConstDependents();
}

Dependency::Dependency(
  ConstParameterEntryList dependees, 
  RCP<ParameterEntry> dependent):
  dependees_(dependees), 
  dependents_(ParameterEntryList(&dependent, &dependent+1))
{
  checkDependeesAndDependents();
  createConstDependents();
}


Dependency::Dependency(
  RCP<const ParameterEntry> dependee, 
  ParameterEntryList dependents):
  dependees_(ConstParameterEntryList(&dependee, &dependee+1)),
  dependents_(dependents)
{
  checkDependeesAndDependents();
  createConstDependents();
}
  
Dependency::Dependency(
  RCP<const ParameterEntry> dependee, 
  RCP<ParameterEntry> dependent):
  dependees_(ConstParameterEntryList(&dependee, &dependee+1)),
  dependents_(ParameterEntryList(&dependent, &dependent+1))
{
  checkDependeesAndDependents();
  createConstDependents();
}


void Dependency::createConstDependents(){
  for(
    ParameterEntryList::iterator it = dependents_.begin();
    it != dependents_.end();
    ++it)
  {
    constDependents_.insert(it->getConst());
  }
}

void Dependency::print(std::ostream& out) const{
  out << "Type: " << getTypeAttributeValue() << std::endl;
  out << "Number of dependees: " << dependees_.size() << std::endl;
  out << "Number of dependents: " << dependents_.size() << std::endl;

}

void Dependency::checkDependeesAndDependents(){
  ConstParameterEntryList::iterator it1 = dependees_.begin(); 
  for(; it1 != dependees_.end(); ++it1){
    TEST_FOR_EXCEPTION((*it1).is_null(),
      InvalidDependencyException,
      "Cannot have a null dependee!" << std::endl << std::endl);
   }

  ParameterEntryList::iterator it2 = dependents_.begin(); 
  for(; it2 != dependents_.end(); ++it2){
    TEST_FOR_EXCEPTION((*it2).is_null(),
      InvalidDependencyException,
      "Cannot have a null dependent!" << std::endl << std::endl);
  }
}

} //namespace Teuchos

