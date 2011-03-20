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

#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_DependencyXMLConverterDB.hpp"


namespace Teuchos {


XMLParameterListReader::XMLParameterListReader()
{;}

RCP<ParameterList> XMLParameterListReader::toParameterList(
  const XMLObject& xml, RCP<DependencySheet> depSheet) const 
{
  TEST_FOR_EXCEPTION(
    xml.getTag() 
    != 
    XMLParameterListWriter::getParameterListTagName(), 
    BadXMLParameterListRootElementException,
    "XMLParameterListReader expected tag " << 
    XMLParameterListWriter::getParameterListTagName()
    <<", found " << xml.getTag());
  RCP<ParameterList> rtn = rcp(new ParameterList);
  IDtoValidatorMap validatorIDsMap;
  int validatorsIndex = 
    xml.findFirstChild(XMLParameterListWriter::getValidatorsTagName());
  if(validatorsIndex != -1){
    convertValidators(xml.getChild(validatorsIndex), validatorIDsMap);
  }
  EntryIDsMap entryIDsMap; 
  convertParameterList(xml, rtn, entryIDsMap, validatorIDsMap);
  
  int dependencyIndex = xml.findFirstChild(
    XMLParameterListWriter::getDependenciesTagName());
  if(dependencyIndex != -1){
    convertDependencies(
      depSheet, 
      xml.getChild(dependencyIndex), 
      entryIDsMap,
      validatorIDsMap);
  }
  return rtn;
}


ParameterList 
XMLParameterListReader::toParameterList(const XMLObject& xml) const
{
  TEST_FOR_EXCEPTION(
    xml.getTag() 
    != 
    XMLParameterListWriter::getParameterListTagName(), 
    BadXMLParameterListRootElementException,
    "XMLParameterListReader expected tag " << 
    XMLParameterListWriter::getParameterListTagName()
    <<", found " << xml.getTag());
  RCP<ParameterList> rtn = rcp(new ParameterList);
  IDtoValidatorMap validatorIDsMap;
  int validatorsIndex = 
    xml.findFirstChild(XMLParameterListWriter::getValidatorsTagName());
  if(validatorsIndex != -1){
    convertValidators(xml.getChild(validatorsIndex), validatorIDsMap);
  }
  EntryIDsMap entryIDsMap; 
  convertParameterList(xml, rtn, entryIDsMap, validatorIDsMap);
  ParameterList toReturn = ParameterList(*rtn);
  return toReturn;
}


void XMLParameterListReader::convertValidators(
  const XMLObject& xml, IDtoValidatorMap& validatorIDsMap) const
{
  std::set<const XMLObject*> validatorsWithPrototypes;
  for (int i=0; i<xml.numChildren(); ++i){
    if (xml.getChild(i).hasAttribute(
      ValidatorXMLConverter::getPrototypeIdAttributeName()))
    {
      validatorsWithPrototypes.insert(&xml.getChild(i));
    }
    else{
      RCP<ParameterEntryValidator> insertedValidator = 
        ValidatorXMLConverterDB::convertXML(
          xml.getChild(i), validatorIDsMap);
      ParameterEntryValidator::ValidatorID xmlID = 
        xml.getChild(i).getRequired<ParameterEntryValidator::ValidatorID>(
        ValidatorXMLConverter::getIdAttributeName());
      testForDuplicateValidatorIDs(xmlID, validatorIDsMap);
      validatorIDsMap.insert(IDtoValidatorMap::IDValidatorPair(
       xmlID,
       insertedValidator));
    }
  }

  for (
    std::set<const XMLObject*>::const_iterator it = 
      validatorsWithPrototypes.begin();
    it!=validatorsWithPrototypes.end();
    ++it)
  {
    RCP<ParameterEntryValidator> insertedValidator =
      ValidatorXMLConverterDB::convertXML(*(*it), validatorIDsMap);
    ParameterEntryValidator::ValidatorID xmlID = 
      (*it)->getRequired<ParameterEntryValidator::ValidatorID>(
         ValidatorXMLConverter::getIdAttributeName());
    testForDuplicateValidatorIDs(xmlID, validatorIDsMap);
    validatorIDsMap.insert(IDtoValidatorMap::IDValidatorPair(
      xmlID, insertedValidator));
  }
}

      
void
XMLParameterListReader::convertParameterList(const XMLObject& xml, 
  RCP<ParameterList> parentList,
  EntryIDsMap& entryIDsMap, const IDtoValidatorMap& validatorIDsMap) const
{
  TEST_FOR_EXCEPTION(
    xml.getTag() != XMLParameterListWriter::getParameterListTagName(), 
    BadParameterListElementException,
    "XMLParameterListReader expected tag " << 
    XMLParameterListWriter::getParameterListTagName()
    <<", found the tag "
    << xml.getTag());
        
  if(xml.hasAttribute(XMLParameterListWriter::getNameAttributeName())){
    parentList->setName(
      xml.getRequired(XMLParameterListWriter::getNameAttributeName()));
  }

  for (int i=0; i<xml.numChildren(); i++) {

      XMLObject child = xml.getChild(i);
      
      TEST_FOR_EXCEPTION(
        child.getTag() != XMLParameterListWriter::getParameterListTagName() 
        &&
        child.getTag() != ParameterEntry::getTagName()
        &&
        child.getTag() != XMLParameterListWriter::getValidatorsTagName()
        &&
        child.getTag() != XMLParameterListWriter::getDependenciesTagName(),
        BadParameterListElementException,
        "XMLParameterListReader expected tag "
        << XMLParameterListWriter::getParameterListTagName() << " or "
        << ParameterEntry::getTagName() << ", but found "
        << child.getTag() << " tag.");


      if(
        child.getTag() == XMLParameterListWriter::getParameterListTagName() 
        || 
        child.getTag() == ParameterEntry::getTagName()
        )
      {
        TEST_FOR_EXCEPTION(
          !child.hasAttribute(XMLParameterListWriter::getNameAttributeName()),
          NoNameAttributeExecption,
          "All child nodes of a ParameterList must have a name attribute!" <<
          std::endl << std::endl);

        const std::string& name =
          child.getRequired(XMLParameterListWriter::getNameAttributeName());
        
        if (child.getTag()==XMLParameterListWriter::getParameterListTagName()) {
          RCP<ParameterList> newList = sublist(parentList, name);
            convertParameterList(child, newList, entryIDsMap, validatorIDsMap);
        }
        else if (child.getTag() == ParameterEntry::getTagName()) {
          parentList->setEntry(
            name, ParameterEntryXMLConverterDB::convertXML(child));
          if(child.hasAttribute(ValidatorXMLConverter::getIdAttributeName())){
            IDtoValidatorMap::const_iterator result = validatorIDsMap.find(
              child.getRequired<ParameterEntryValidator::ValidatorID>(
                ValidatorXMLConverter::getIdAttributeName()));
            TEST_FOR_EXCEPTION(result == validatorIDsMap.end(), 
              MissingValidatorDefinitionException,
              "Could not find validator with id: "
              << child.getRequired(
                ValidatorXMLConverter::getIdAttributeName())
              << std::endl << 
              "Bad Parameter: " << name << std::endl << std::endl);
            parentList->getEntryRCP(name)->setValidator(result->second);
        }  
      } 
      if(child.hasAttribute(ParameterEntryXMLConverter::getIdAttributeName())){
        insertEntryIntoMap(child, parentList->getEntryRCP(name), entryIDsMap);
      }
    }  
  }
}

void XMLParameterListReader::testForDuplicateValidatorIDs(
  ParameterEntryValidator::ValidatorID potentialNewID,
  const IDtoValidatorMap& currentMap) const
{
  TEST_FOR_EXCEPTION(currentMap.find(potentialNewID) != currentMap.end(),
  DuplicateValidatorIDsException,
  "Validators with duplicate ids found!" << std::endl <<
  "Bad ID: " << potentialNewID);
}

void XMLParameterListReader::convertDependencies(
  RCP<DependencySheet> depSheet, 
  const XMLObject& xml, 
  const EntryIDsMap& entryIDsMap,
  const IDtoValidatorMap& validatorIDsMap) const
{
  if(xml.hasAttribute(DependencySheet::getNameAttributeName())){
    depSheet->setName(
      xml.getAttribute(DependencySheet::getNameAttributeName()));
  }
  for(int i = 0; i < xml.numChildren(); ++i){
    RCP<Dependency> currentDep = DependencyXMLConverterDB::convertXML(
      xml.getChild(i), 
      entryIDsMap, 
      validatorIDsMap);
    depSheet->addDependency(currentDep);
  }
}

void XMLParameterListReader::insertEntryIntoMap(
  const XMLObject& xmlObj,
  RCP<ParameterEntry> entryToInsert,
  EntryIDsMap& entryIDsMap) const
{
  if(xmlObj.hasAttribute(ParameterEntryXMLConverter::getIdAttributeName()))
  {
     ParameterEntry::ParameterEntryID xmlID = 
       xmlObj.getRequired<ParameterEntry::ParameterEntryID>(
          ParameterEntryXMLConverter::getIdAttributeName());
     TEST_FOR_EXCEPTION(entryIDsMap.find(xmlID) != entryIDsMap.end(),
        DuplicateParameterIDsException,
       "Parameters/ParameterList with duplicate ids found!" << std::endl <<
       "Bad ID: " << xmlID << std::endl << std::endl);
     entryIDsMap.insert(EntryIDsMap::value_type(xmlID, entryToInsert));
  }
}


} // namespace Teuchos

