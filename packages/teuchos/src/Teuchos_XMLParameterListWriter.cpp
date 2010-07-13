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

using namespace Teuchos;

XMLParameterListWriter::XMLParameterListWriter()
{;}

XMLObject XMLParameterListWriter::toXML(const ParameterList& p) const
{
	XMLObject toReturn(getParameterListAspectsTagName());
	WriterValidatorIDMap validatorIDMap;
	int validatorIDCounter = 0;
	XMLObject parameterLists(getParameterListsTagName());
	parameterLists.addChild(convertParameterList(p, validatorIDMap, validatorIDCounter));
	WriterValidatorIDMap::const_iterator it = validatorIDMap.begin();
	XMLObject validators(getValidatorsTagName());
	for(;it != validatorIDMap.end(); ++it){
		XMLObject currentValidator = ValidatorXMLConverterDB::getConverter(*(it->first))->fromValidatortoXML(it->first)
		currentValidator.addInt(getValidatorIdAttributeName(), it->second);
		validators.addChild(currentValidator);
	}
	toReturn.addChild(validators);
	toReturn.addChild(parameterLists);
	return toReturn;
}

XMLObject XMLParameterListWriter::convertParameterList(
	const ParameterList& p,
	WriterValidatorIDMap& validatorIDMap,
	int& validatorIDCounter) const
{
  XMLObject rtn(getParameterListTagName());
  rtn.addAttribute(getNameAttributeName(), p.name());
  
  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i)
    {
     const ParameterEntry& entry = p.entry(i);
     if (entry.isList())
       {
          rtn.addChild(convertParameterList(getValue<ParameterList>(entry), validatorIDMap, validatorIDCounter));
       }
     else
	   {
          const std::string& name = p.name(i);
		  RCP<const ParameterEntryXMLConverter> converter = ParameterEntryXMLConverterDB::getConverter(entry);
		  XMLObject parameterEntryXML = converter->fromParameterEntrytoXML(entry, name);
		  if(!entry.validator().is_null()){
			WriterValidatorIDMap::iterator result = validatorIDMap.find(entry.validator().get());
			if(result == validatorIDMap.end()){
			  validatorIDMap.insert(WriterValidatorIDPair(entry.validator().get(), validatorIDCounter));
			  parameterEntryXML.addInt(getValidatorIdAttributeName(), validatorIDCounter);
			  ++validatorIDCounter;
			}
			else{
			  parameterEntryXML.addInt(getValidatorIdAttributeName(), result->second);	
			}
		  }
		  rtn.addChild(parameterEntryXML);
	   }
    }
  return rtn;
}

