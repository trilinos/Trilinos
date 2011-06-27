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

