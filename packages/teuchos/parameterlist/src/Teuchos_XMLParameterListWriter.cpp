// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  const ParameterList& /* p */, ValidatortoIDMap& validatorIDsMap) const
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

