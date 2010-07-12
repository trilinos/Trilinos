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


using namespace Teuchos;

XMLParameterListReader::XMLParameterListReader()
{;}

ParameterList XMLParameterListReader::toParameterList(const XMLObject& xml) const{
	TEST_FOR_EXCEPTION(xml.getTag() != XMLParameterListWriter::getParameterListAspectsTagName(), std::runtime_error,
		"XMLParameterListReader expected tag " << XMLParameterListWriter::getParameterListAspectsTagName() <<", found "
		<< xml.getTag());
	ReaderValidatorIDMap validators;
	ParameterList rtn;
	for(int i =0; i<xml.numChildren(); ++i){
		if(xml.getChild(i).getTag() == XMLParameterListWriter::getValidatorsTagName()){
			convertValidators(xml.getChild(i), validators);
		}
	}
	for(int i =0; i<xml.numChildren(); ++i){
		if(xml.getChild(i).getTag() == XMLParameterListWriter::getParameterListsTagName()){
			rtn = convertParameterList(xml.getChild(i).getChild(0), validators);
		}
	}
	return rtn;
}

void XMLParameterListReader::convertValidators(const XMLObject& xml, ReaderValidatorIDMap& validators) const
{
	for(int i=0; i<xml.numChildren(); ++i){
		int currentID = xml.getChild(i).getRequiredInt(XMLParameterListWriter::getValidatorIdAttributeName());
		RCP<ParameterEntryValidator> currentValidator = ValidatorXMLConverterDB::getConverter(xml.getChild(i))->fromXMLtoValidator(xml.getChild(i));
		validators.insert(ReaderValidatorIDPair(currentID, currentValidator));
	}
}
			
ParameterList XMLParameterListReader::convertParameterList(const XMLObject& xml, const ReaderValidatorIDMap& validators) const
{
  TEST_FOR_EXCEPTION(xml.getTag() != XMLParameterListWriter::getParameterListTagName(), std::runtime_error,
                     "XMLParameterListReader expected tag " << XMLParameterListWriter::getParameterListTagName() <<", found "
                     << xml.getTag());

  ParameterList rtn;

  for (int i=0; i<xml.numChildren(); i++)
    {
      XMLObject child = xml.getChild(i);

      TEST_FOR_EXCEPTION( (child.getTag() != XMLParameterListWriter::getParameterListTagName() 
                           && child.getTag() != ParameterEntry::getTagName()), 
                         std::runtime_error,
                         "XMLParameterListReader expected tag "
                         << XMLParameterListWriter::getParameterListsTagName() << " or "
						 << ParameterEntry::getTagName() << ", found "
                         << child.getTag());

      if (child.getTag()==XMLParameterListWriter::getParameterListTagName())
        {
          TEST_FOR_EXCEPTION( !child.hasAttribute(XMLParameterListWriter::getNameAttributeName()), 
                         std::runtime_error,
                         XMLParameterListWriter::getParameterListTagName() <<" tags must "
						 "have a " << XMLParameterListWriter::getNameAttributeName() << " attribute"
                         << child.getTag());

		  const std::string& name = child.getRequired(XMLParameterListWriter::getNameAttributeName());

          ParameterList sublist = convertParameterList(child, validators);
          sublist.setName(name);

          rtn.set(name, sublist);
        }
      else
        {
          const std::string& name = child.getRequired(XMLParameterListWriter::getNameAttributeName());
			RCP<const ParameterEntryXMLConverter> converter = ParameterEntryXMLConverterDB::getConverter(child);
			ParameterEntry parameter = converter->fromXMLtoParameterEntry(child);
			if(child.hasAttribute(XMLParameterListWriter::getValidatorIdAttributeName())){
				ReaderValidatorIDMap::const_iterator result = validators.find(child.getRequiredInt(XMLParameterListWriter::getValidatorIdAttributeName()));
				if(result != validators.end()){
					parameter.setValidator(result->second);
				}
				else{
					throw std::runtime_error("Could not find validator with id: " + toString(child.getRequiredInt(XMLParameterListWriter::getValidatorIdAttributeName())));
				}
			}	
			rtn.setEntry(name, converter->fromXMLtoParameterEntry(child));
        } 
    }
  return rtn;
}


