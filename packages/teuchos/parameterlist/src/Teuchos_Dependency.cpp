// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    TEUCHOS_TEST_FOR_EXCEPTION((*it1).is_null(),
      InvalidDependencyException,
      "Cannot have a null dependee!" << std::endl << std::endl);
   }

  ParameterEntryList::iterator it2 = dependents_.begin();
  for(; it2 != dependents_.end(); ++it2){
    TEUCHOS_TEST_FOR_EXCEPTION((*it2).is_null(),
      InvalidDependencyException,
      "Cannot have a null dependent!" << std::endl << std::endl);
  }
}

} //namespace Teuchos

