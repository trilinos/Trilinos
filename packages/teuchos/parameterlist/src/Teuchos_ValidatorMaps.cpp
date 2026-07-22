// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ValidatorMaps.hpp"

namespace Teuchos{


void IDtoValidatorMap::insert(IDValidatorPair toInsert)
{
  validatorMap.insert(toInsert);
}


IDtoValidatorMap::const_iterator IDtoValidatorMap::find(int id) const
{
  return validatorMap.find(id);
}


IDtoValidatorMap::const_iterator IDtoValidatorMap::begin() const
{
  return validatorMap.begin();
}


IDtoValidatorMap::const_iterator IDtoValidatorMap::end() const
{
  return validatorMap.end();
}


ValidatortoIDMap::ValidatortoIDMap():counter(0)
{}

void ValidatortoIDMap::insert(RCP<const ParameterEntryValidator> toInsert)
{
  const_iterator result = validatorMap.find(toInsert);
  if(result == validatorMap.end()){
    validatorMap.insert(ValidatorIDPair(toInsert, counter));
    ++counter;
  }
}


ValidatortoIDMap::const_iterator ValidatortoIDMap::find(
  const RCP<const ParameterEntryValidator> validator) const
{
  return validatorMap.find(validator);
}


ValidatortoIDMap::const_iterator ValidatortoIDMap::begin() const
{
  return validatorMap.begin();
}


ValidatortoIDMap::const_iterator ValidatortoIDMap::end() const
{
  return validatorMap.end();
}


} // namespace Teuchos
