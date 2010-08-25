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

#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"


namespace Teuchos{


XMLParameterListWriter::XMLParameterListWriter()
{;}


XMLObject XMLParameterListWriter::toXML(const ParameterList& p) const
{
  XMLObject validators = convertValidators(p);
  XMLObject toReturn = convertParameterList(p);
  toReturn.addChild(validators);
  return toReturn;
}

void XMLParameterListWriter::buildValidatorSet(
  const ParameterList& p,
  ValidatorSet& validatorSet) const
{
  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i) {
    const ParameterEntry& entry = p.entry(i);
    if(entry.isList()){
      buildValidatorSet(
        getValue<ParameterList>(entry),
        validatorSet);
    }
    else if(nonnull(entry.validator())){
      validatorSet.insert(entry.validator());
    }
  }
}


XMLObject XMLParameterListWriter::convertValidators(
  const ParameterList& p) const
{
  XMLObject validators(getValidatorsTagName());
  ValidatorSet validatorSet;
  buildValidatorSet(p, validatorSet);
  for(
    ValidatorSet::const_iterator it = validatorSet.begin();
    it != validatorSet.end();
    ++it)
  {
    validators.addChild(
      ValidatorXMLConverterDB::convertValidator(*it, validatorSet));
  }
  return validators;
}


XMLObject XMLParameterListWriter::convertParameterList(
  const ParameterList& p) const
{
  XMLObject rtn(getParameterListTagName());
  rtn.addAttribute(getNameAttributeName(), p.name());
  
  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i){
    RCP<const ParameterEntry> entry = p.getEntryRCP(i->first);
    if(entry->isList()){
      rtn.addChild(convertParameterList(getValue<ParameterList>(entry)));
    }
    else{
      rtn.addChild(ParameterEntryXMLConverterDB::convertEntry(
        entry, p.name(i)));
    }
  }
  return rtn;
}


} // namespace Teuchos
