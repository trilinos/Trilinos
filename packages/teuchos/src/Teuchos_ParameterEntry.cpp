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


#include "Teuchos_ParameterEntry.hpp" // class definition
#include "Teuchos_ParameterList.hpp"	 // for sublists


namespace Teuchos {


ParameterEntry::ParameterEntry() : 
  isUsed_(false),
  isDefault_(false)
{
  addParameterEntryToIDMap(this);
}


ParameterEntry::ParameterEntry(const ParameterEntry& source)
{
  operator=(source);
}


ParameterEntry& ParameterEntry::operator=(const ParameterEntry& source)
{
  if (&source == this)
    return *this;

  val_ = source.val_;
  isUsed_ = source.isUsed_;
  isDefault_ = source.isDefault_;
  docString_ = source.docString_;
  validator_ = source.validator_;

  addParameterEntryToIDMap(this);

  return *this;
}


void ParameterEntry::setAnyValue(
  const any &value_in, bool isDefault_in
  )
{
  val_ = value_in;
  isDefault_ = isDefault_in;
  validator_ = null;
  isUsed_ = false;
  docString_ = "";
}


void ParameterEntry::setValidator(
  RCP<const ParameterEntryValidator> const& validator_in
  )
{
  validator_ = validator_in;
}


void ParameterEntry::setDocString(const std::string &docString_in)
{
  docString_ = docString_in;
}


ParameterList& ParameterEntry::setList(
  bool isDefault_in, const std::string &docString_in
  )
{
  val_ = ParameterList();
  isDefault_ = isDefault_in;
  isUsed_ = true;
  docString_ = docString_in;
  return any_cast<ParameterList>( val_ );
}


bool ParameterEntry::isList() const
{
  return ( val_.empty() ? false : val_.type() == typeid(ParameterList) );
}


std::ostream& ParameterEntry::leftshift(std::ostream& os, bool printFlags) const
{
  if( !this->isList() ) os << val_;

  if(printFlags) {
    if (isDefault_)
      os << "   [default]";
    else if (!isUsed_)
      os << "   [unused]";
  }

  return os;
}

RCP<ParameterEntry> ParameterEntry::getParameterEntry(ParameterEntryID id){
  IDToParameterEntryMap::iterator result = find(
    masterIDMap.begin(), masterIDMap.end(), id);
  result != masterIDMap.end() ? return result->second : return null;
}

ParameterEntryID ParameterEntry::getParameterEntryID(RCP<ParameterEntry> entry){
  ParameterEntryToIDMap::iterator result = find(
    masterParameterEntryMap.begin(), masterParameterEntryMap.end(), entry);
  result != masterParameterEntryMap.end() ? 
    return it->second : return OrdinalTraits<ParameterEntryID>::invalid();
}

// private


void ParameterEntry::reset()
{
  //delete val_;
  isUsed_ = false;
  isDefault_ = false;
}

void ParameterEntry::incrementMasterCounter(){
  masterIDCounter++;
  if(masterIDMap.find(masterIDCounter) != masterIDMap.end()){
    incrementMasterCounter;
  }
}

ParameterEntryID ParameterEntry::getAvailableID(){
  ParameterEntryID toReturn;
  if(masterFreeIDs.size() != 0)
    toReturn = masterFreeIDs.back();
    masterFreeIDs.pop_back();
    return toReturn;
  }
  else{
    toReturn = masterIDCounter;
    incrementMasterCounter();
  }
  return toReturn;
}

void ParameterEntry::addParameterEntryToIDMap(ParameterEntry* entryToAdd){
  RCP<ParameterEntry> toInsert = rcp(entryToAdd, false);
  TEST_FOR_EXCEPTION(
  masterParameterEntryMap.find(toInsert)
  !=
  masterParameterEntryMap.end(),
  DuplicateParameterEntryException,
  "Error: parameter entry has already been added to the parameter entry maps!" <<
  std::endl << std::endl);

  ParameterEntryID insertionID = getAvailableID();

  masterIDMap.insert(value_type(insertionID, toInsert));
  masterParameterEntryMap.insert(value_type(toInsert, insertionID));
}

void ParameterEntry::addParameterEntryToIDMap(
  ParameterEntry* entry, ParameterEntryID idToUse)
{
  TEST_FOR_EXCEPTION(
  masterIDMap.find(idToUse)
  !=
  masterIDMap.end(),
  DuplicateParameterEntryIDException,
  "Error: can't create a ParameterEntry with the ID" <<
  insertionID << ". That ID is already being used to track another " <<
  "ParameterEntry!" << std::endl << std::endl);
  masterParameterEntryMap.insert( value_type(rcp(entry, false), idToUse));
  masterIDMap.insert(value_type(idToUse, rcp(entry, false)));
}

void ParameterEntry::removeParameterEntryFromMasterMaps(Parameter* entryToRemove){
  ParameterEntryToIDMap::iterator toRemove = 
    masterParameterEntryMap.find(rcp(entryToRemove, false));
  ParameterEntryID idToFree = toRemove->second;
  masterParameterEntryMap.remove(toRemove);
  masterIDMap.erase(idToFree);
  masterFreeIDs.push_back(idToFree);
}


} // namespace Teuchos

