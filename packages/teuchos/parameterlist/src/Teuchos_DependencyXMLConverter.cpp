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

#include "Teuchos_DependencyXMLConverter.hpp"
#include "Teuchos_XMLDependencyExceptions.hpp"


namespace Teuchos{


RCP<Dependency>
DependencyXMLConverter::fromXMLtoDependency(
  const XMLObject& xmlObj,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap,
  const IDtoValidatorMap& validatorIDsMap) const
{
  Dependency::ConstParameterEntryList dependees;
  Dependency::ParameterEntryList dependents;

  TEUCHOS_TEST_FOR_EXCEPTION(xmlObj.findFirstChild(getDependeeTagName()) == -1,
    MissingDependeesException,
    "Could not find any dependees for a dependency!" 
    <<std::endl <<std::endl);

  TEUCHOS_TEST_FOR_EXCEPTION(xmlObj.findFirstChild(getDependentTagName()) == -1,
    MissingDependentsException,
    "Could not find any dependents for a dependency!" 
    <<std::endl <<std::endl);

  for(int i=0; i<xmlObj.numChildren(); ++i){
    XMLObject child = xmlObj.getChild(i);
    if(child.getTag() == getDependeeTagName()){
      ParameterEntry::ParameterEntryID dependeeID = 
        child.getRequired<ParameterEntry::ParameterEntryID>(
          getParameterIdAttributeName());

      TEUCHOS_TEST_FOR_EXCEPTION(entryIDsMap.find(dependeeID) == entryIDsMap.end(),
        MissingDependeeException,
        "Can't find a Dependee ParameterEntry associated with the ID: " <<
        dependeeID << std::endl << std::endl);
      dependees.insert(entryIDsMap.find(dependeeID)->second);
    }
    else if(child.getTag() == getDependentTagName()){
      ParameterEntry::ParameterEntryID dependentID = 
        child.getRequired<ParameterEntry::ParameterEntryID>(
          getParameterIdAttributeName());

      TEUCHOS_TEST_FOR_EXCEPTION(entryIDsMap.find(dependentID) == entryIDsMap.end(),
        MissingDependentException,
        "Can't find a Dependent ParameterEntry associated with the ID: " <<
        dependentID << std::endl << std::endl);
      dependents.insert(entryIDsMap.find(dependentID)->second);
    }
  }

  return convertXML(
    xmlObj, dependees, dependents, entryIDsMap, validatorIDsMap);
}

XMLObject
DependencyXMLConverter::fromDependencytoXML(
  const RCP<const Dependency> dependency,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
  ValidatortoIDMap& validatorIDsMap) const
{
  XMLObject toReturn(Dependency::getXMLTagName());

  toReturn.addAttribute(getTypeAttributeName(), 
    dependency->getTypeAttributeValue());

  Dependency::ConstParameterEntryList::const_iterator it =
    dependency->getDependees().begin();

  for(;it != dependency->getDependees().end(); ++it){
    XMLObject currentDependee(getDependeeTagName());
    TEUCHOS_TEST_FOR_EXCEPTION(entryIDsMap.find(*it) == entryIDsMap.end(),
      MissingDependeeException,
      "Can't find the Dependee of a dependency in the given " <<
      "EntryIDsMap. Occured when converting " <<
      "to XML" << std::endl << std::endl);
    currentDependee.addAttribute<ParameterEntry::ParameterEntryID>(
      getParameterIdAttributeName(), entryIDsMap.find(*it)->second);
    toReturn.addChild(currentDependee);
  }

  it = dependency->getDependents().begin();
  for(; it != dependency->getDependents().end(); ++it){
    XMLObject currentDependent(getDependentTagName());
    TEUCHOS_TEST_FOR_EXCEPTION(entryIDsMap.find(*it) == entryIDsMap.end(),
      MissingDependentException,
      "Can't find the Dependent of a dependency in the given " <<
      "ValidatordIDsMap.. Occured when converting " <<
      "to XML" << std::endl << std::endl);
    currentDependent.addAttribute<ParameterEntry::ParameterEntryID>(
      getParameterIdAttributeName(), entryIDsMap.find(*it)->second);
    toReturn.addChild(currentDependent);
  }

  convertDependency(dependency, toReturn, entryIDsMap, validatorIDsMap);

  return toReturn;
}


} // namespace Teuchos

