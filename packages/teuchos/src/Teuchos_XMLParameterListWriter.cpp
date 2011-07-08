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

#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_DependencyXMLConverterDB.hpp"


namespace Teuchos{


XMLParameterListWriter::XMLParameterListWriter()
{;}


XMLObject 
XMLParameterListWriter::toXML(
  const ParameterList& p, 
  RCP<const DependencySheet> depSheet) const
{
  EntryIDsMap entryIDsMap;
  ValidatortoIDMap validatorIDsMap;
  ParameterEntry::ParameterEntryID peIDCounter = 0;

  //We build an initial map full of validators that are located in the 
  //parameter list. That way we can convert the parameter entries.
  buildInitialValidatorMap(p, validatorIDsMap);

  XMLObject toReturn = 
    convertParameterList(p, peIDCounter, entryIDsMap, validatorIDsMap);
  toReturn.addAttribute(getNameAttributeName(), p.name());

  if(!depSheet.is_null()){
    XMLObject deps = 
      convertDependencies(depSheet, entryIDsMap, validatorIDsMap);
    toReturn.addChild(deps);
  }

  //Validators must be done after depencneies because dependencies might add
  //entries to the validator map. KLN 09/20/2010
  XMLObject validators = convertValidators(p, validatorIDsMap);
  toReturn.addChild(validators);

  return toReturn;
}

void XMLParameterListWriter::buildInitialValidatorMap(
  const ParameterList& p,
  ValidatortoIDMap& validatorIDsMap) const
{
  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i) {
    const ParameterEntry& entry = p.entry(i);
    if(entry.isList()){
      buildInitialValidatorMap(
        getValue<ParameterList>(entry),
        validatorIDsMap);
    }
    else if(nonnull(entry.validator())){
      validatorIDsMap.insert(entry.validator());
    }
  }
}


XMLObject XMLParameterListWriter::convertValidators(
  const ParameterList& p, ValidatortoIDMap& validatorIDsMap) const
{
  XMLObject validators(getValidatorsTagName());
  for(
    ValidatortoIDMap::const_iterator it = validatorIDsMap.begin();
    it != validatorIDsMap.end();
    ++it)
  {
    validators.addChild(
      ValidatorXMLConverterDB::convertValidator(it->first, validatorIDsMap));
  }
  return validators;
}


XMLObject XMLParameterListWriter::convertParameterList(
  const ParameterList& p,
  ParameterEntry::ParameterEntryID& idCounter,
  EntryIDsMap& entryIDsMap,
  const ValidatortoIDMap& validatorIDsMap) const
{
  XMLObject rtn(getParameterListTagName());
  
  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i){
    RCP<const ParameterEntry> entry = p.getEntryRCP(i->first);
    if(entry->isList()){
      XMLObject newPL = convertParameterList(
        getValue<ParameterList>(entry), 
        idCounter, 
        entryIDsMap,
        validatorIDsMap);
      newPL.addAttribute(
        getNameAttributeName(), p.name(i));
      newPL.addAttribute(
        ParameterEntryXMLConverter::getIdAttributeName(), idCounter);
      entryIDsMap.insert(EntryIDsMap::value_type(entry, idCounter));
      rtn.addChild(newPL);
      ++idCounter;
    }
    else{
      rtn.addChild(ParameterEntryXMLConverterDB::convertEntry(
        entry, p.name(i), idCounter, validatorIDsMap));
      entryIDsMap.insert(EntryIDsMap::value_type(entry, idCounter));
      ++idCounter;
    }
  }
  return rtn;
}

XMLObject 
XMLParameterListWriter::convertDependencies(
  RCP<const DependencySheet> depSheet,
  const EntryIDsMap& entryIDsMap,
  ValidatortoIDMap& validatorIDsMap) const
{
  XMLObject toReturn(getDependenciesTagName());
  toReturn.addAttribute(
    DependencySheet::getNameAttributeName(), 
    depSheet->getName()
  );

  for(
    DependencySheet::DepSet::const_iterator it = depSheet->depBegin();
    it != depSheet->depEnd();
    ++it)
  {
    toReturn.addChild(DependencyXMLConverterDB::convertDependency(
      *it, entryIDsMap, validatorIDsMap));
  }
  return toReturn;
}


} // namespace Teuchos

