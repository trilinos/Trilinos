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
#include "Teuchos_OrdinalTraits.hpp"


namespace Teuchos {


ParameterEntry::ParameterEntry() : 
  isUsed_(false),
  isDefault_(false)
{
  addParameterEntryToMasterMaps(this);
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

  addParameterEntryToMasterMaps(this);

  return *this;
}

ParameterEntry::~ParameterEntry()
{
  removeParameterEntryFromMasterMaps(this);
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
  IDToParameterEntryMap::iterator result = getMasterIDMap().find(id);
  return result != getMasterIDMap().end() ? result->second : null;
}

ParameterEntry::ParameterEntryID 
ParameterEntry::getParameterEntryID(RCP<const ParameterEntry> entry){
  ParameterEntryToIDMap::iterator result =
    getMasterParameterEntryMap().find(entry);
  return result != getMasterParameterEntryMap().end() ? 
    result->second : OrdinalTraits<ParameterEntryID>::invalid();
}

// private


void ParameterEntry::reset()
{
  //delete val_;
  isUsed_ = false;
  isDefault_ = false;
}

void ParameterEntry::incrementMasterCounter(){
  getMasterIDCounter()++;
  if(getMasterIDMap().find(getMasterIDCounter()) != getMasterIDMap().end()){
    incrementMasterCounter();
  }
}

ParameterEntry::ParameterEntryID ParameterEntry::getAvailableID(){
  ParameterEntryID toReturn;
  if(getMasterFreeIDs().size() != 0){
    toReturn = getMasterFreeIDs().back();
    getMasterFreeIDs().pop_back();
    return toReturn;
  }
  else{
    toReturn = getMasterIDCounter();
    incrementMasterCounter();
  }
  return toReturn;
}

void ParameterEntry::addParameterEntryToMasterMaps(
  ParameterEntry* entryToAdd)
{
  RCP<ParameterEntry> toInsert = rcpFromUndefRef(*entryToAdd);

  ParameterEntryID insertionID = getAvailableID();

  getMasterIDMap().insert(
    IDParameterEntryPair(insertionID, toInsert));
  getMasterParameterEntryMap().insert(
    ParameterEntryIDPair(toInsert.getConst(), insertionID));
}

void ParameterEntry::addParameterEntryToMasterMaps(
  ParameterEntry* entry, ParameterEntryID idToUse)
{
  TEST_FOR_EXCEPTION(
  getMasterIDMap().find(idToUse)
  !=
  getMasterIDMap().end(),
  DuplicateParameterEntryIDException,
  "Error: can't create a ParameterEntry with the ID" <<
  idToUse << ". That ID is already being used to track another " <<
  "ParameterEntry!" << std::endl << std::endl);
  getMasterParameterEntryMap().insert(
    ParameterEntryIDPair(rcpFromUndefRef(*entry).getConst(), idToUse));
  getMasterIDMap().insert(IDParameterEntryPair(idToUse, rcpFromUndefRef(*entry)));
}

void ParameterEntry::removeParameterEntryFromMasterMaps(
  ParameterEntry* entryToRemove)
{
  ParameterEntryToIDMap::iterator toRemove = 
    getMasterParameterEntryMap().find(rcpFromUndefRef(*entryToRemove));
  ParameterEntryID idToFree = toRemove->second;
  getMasterParameterEntryMap().erase(toRemove->first);
  getMasterIDMap().erase(idToFree);
  getMasterFreeIDs().push_back(idToFree);
}


} // namespace Teuchos


