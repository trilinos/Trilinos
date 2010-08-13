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


namespace Teuchos {


XMLParameterListReader::XMLParameterListReader()
{;}


ParameterList XMLParameterListReader::toParameterList(const XMLObject& xml) const
{
  TEST_FOR_EXCEPTION(xml.getTag() != XMLParameterListWriter::getParameterListTagName(), 
    BadXMLParameterListRootElementException,
    "XMLParameterListReader expected tag " << XMLParameterListWriter::getParameterListTagName()
    <<", found " << xml.getTag());
  ParameterList rtn;
  for(int i =0; i<xml.numChildren(); ++i){
    if(xml.getChild(i).getTag() == XMLParameterListWriter::getValidatorsTagName()){
      convertValidators(xml.getChild(i));
    }
  }

  rtn = convertParameterList(xml);
  return rtn;
}


void XMLParameterListReader::convertValidators(const XMLObject& xml) const
{
  std::set<const XMLObject*> validatorsWithPrototypes;
  for (int i=0; i<xml.numChildren(); ++i){
    if (xml.getChild(i).hasAttribute(
      ValidatorXMLConverter::getPrototypeIdAttributeName()))
    {
      validatorsWithPrototypes.insert(&xml.getChild(i));
    }
    else{
      ValidatorXMLConverterDB::convertXML(xml.getChild(i));
    }
  }
  for (
    std::set<const XMLObject*>::const_iterator it = validatorsWithPrototypes.begin();
    it!=validatorsWithPrototypes.end();
    ++it)
  {
    ValidatorXMLConverterDB::convertXML(*(*it));
  }
}

      
ParameterList XMLParameterListReader::convertParameterList(
  const XMLObject& xml) const
{
  TEST_FOR_EXCEPTION(xml.getTag() != XMLParameterListWriter::getParameterListTagName(), 
    BadParameterListElementException,
    "XMLParameterListReader expected tag " << XMLParameterListWriter::getParameterListTagName()
    <<", found the tag "
    << xml.getTag());

  ParameterList rtn;

  for (int i=0; i<xml.numChildren(); i++) {

      XMLObject child = xml.getChild(i);
      
      TEST_FOR_EXCEPTION(
        (
          child.getTag() != XMLParameterListWriter::getParameterListTagName() 
          &&
          child.getTag() != ParameterEntry::getTagName()
          )
        &&
        child.getTag() != XMLParameterListWriter::getValidatorsTagName(),
        BadParameterListElementException,
        "XMLParameterListReader expected tag "
        << XMLParameterListWriter::getParameterListTagName() << " or "
        << ParameterEntry::getTagName() << ", but found "
        << child.getTag() << " tag.");
      
      if (child.getTag()==XMLParameterListWriter::getParameterListTagName()) {
        TEST_FOR_EXCEPTION(
          !child.hasAttribute(XMLParameterListWriter::getNameAttributeName()), 
          NoNameAttributeExecption,
          XMLParameterListWriter::getParameterListTagName() <<" tags must "
          "have a " << XMLParameterListWriter::getNameAttributeName() << " attribute"
          << child.getTag());
        
        const std::string& name =
          child.getRequired(XMLParameterListWriter::getNameAttributeName());

        ParameterList sublist = convertParameterList(child);
        sublist.setName(name);

        rtn.set(name, sublist);
      }
      else if (child.getTag() == ParameterEntry::getTagName()) {

        TEST_FOR_EXCEPTION(
          !child.hasAttribute(XMLParameterListWriter::getNameAttributeName()), 
          NoNameAttributeExecption,
          ParameterEntry::getTagName() <<" tags must " <<
          "have a " << 
          XMLParameterListWriter::getNameAttributeName() << 
          " attribute" << std::endl << "Bad Parameter: " << 
          child.getAttribute(XMLParameterListWriter::getNameAttributeName()) <<
          std::endl << std::endl);

        const std::string& name =
          child.getRequired(XMLParameterListWriter::getNameAttributeName());
        ParameterEntry parameter = 
          ParameterEntryXMLConverterDB::convertXML(child);
        if(child.hasAttribute(ValidatorXMLConverter::getIdAttributeName())){
          RCP<ParameterEntryValidator> result =
            ParameterEntryValidator::getValidator(
              child.getRequiredInt(ValidatorXMLConverter::getIdAttributeName()));
          TEST_FOR_EXCEPTION(result.is_null(), 
            MissingValidatorDefinitionException,
            "Could not find validator with id: "
            << child.getRequiredInt(ValidatorXMLConverter::getIdAttributeName())
            << std::endl << 
            "Bad Parameter: " << name << std::endl << std::endl);
          parameter.setValidator(result);
        }  
        rtn.setEntry(name, parameter);
     } 
  }
  return rtn;
}


} // namespace Teuchos

