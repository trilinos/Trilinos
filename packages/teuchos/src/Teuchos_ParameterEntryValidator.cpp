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


#include "Teuchos_ParameterEntryValidator.hpp" // class definition
#include "Teuchos_OrdinalTraits.hpp"


namespace Teuchos {

ParameterEntryValidator::ParameterEntryValidator(){
  addValidatorToMasterMaps(this);
}

ParameterEntryValidator::ParameterEntryValidator(ValidatorID validatorID){
  if(validatorID != OrdinalTraits<ValidatorID>::invalid()){
    addValidatorToMasterMaps(this, validatorID);
  }
}

ParameterEntryValidator::~ParameterEntryValidator(){
  removeValidatorFromMasterMaps(this);
}

RCP<ParameterEntryValidator> 
ParameterEntryValidator::getValidator(ValidatorID id)
{
  IDToValidatorMap::iterator result = getMasterIDMap().find(id);
  return result != getMasterIDMap().end() ? rcpFromRef(*(result->second)) : null;
}

ParameterEntryValidator::ValidatorID
ParameterEntryValidator::getValidatorID(
  RCP<const ParameterEntryValidator> validator)
{
  ValidatorToIDMap::iterator result =
    getMasterValidatorMap().find(validator.get());
  return result != getMasterValidatorMap().end() ? 
    result->second : OrdinalTraits<ValidatorID>::invalid();
}

void ParameterEntryValidator::printKnownValidators(std::ostream &out){
  for(
    IDToValidatorMap::iterator it = getMasterIDMap().begin();
    it != getMasterIDMap().end();
    ++it)
  {
    out << it->first; 
  }
}

//private 

void ParameterEntryValidator::incrementMasterCounter(){
  getMasterIDCounter()++;
  if(getMasterIDMap().find(getMasterIDCounter()) != getMasterIDMap().end()){
    incrementMasterCounter();
  }
}

ParameterEntryValidator::ValidatorID ParameterEntryValidator::getAvailableID(){
  ValidatorID toReturn;
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

void ParameterEntryValidator::addValidatorToMasterMaps(
  ParameterEntryValidator* validatorToAdd)
{

  ValidatorID insertionID = getAvailableID();

  getMasterIDMap().insert(
    IDValidatorPair(insertionID, validatorToAdd));
  getMasterValidatorMap().insert(
    ValidatorIDPair(validatorToAdd, insertionID));
}

void ParameterEntryValidator::addValidatorToMasterMaps(
  ParameterEntryValidator* entry, ValidatorID idToUse)
{
  TEST_FOR_EXCEPTION(
  getMasterIDMap().find(idToUse)
  !=
  getMasterIDMap().end(),
  DuplicateValidatorIDException,
  "Error: can't create a ParameterEntryValidator with the ID" <<
  idToUse << ". That ID is already being used to track another " <<
  "ParameterEntryValidator!" << std::endl << std::endl);
  getMasterValidatorMap().insert(
    ValidatorIDPair(entry, idToUse));
  getMasterIDMap().insert(IDValidatorPair(idToUse, entry));
}

void ParameterEntryValidator::removeValidatorFromMasterMaps(
  ParameterEntryValidator* validatorToRemove)
{
  ValidatorToIDMap::iterator toRemove = 
    getMasterValidatorMap().find(validatorToRemove);
  ValidatorID idToFree = toRemove->second;
  getMasterValidatorMap().erase(toRemove->first);
  getMasterIDMap().erase(idToFree);
  getMasterFreeIDs().push_back(idToFree);
}


} // namespace Teuchos

