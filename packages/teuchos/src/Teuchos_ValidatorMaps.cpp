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
