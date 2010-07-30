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

#include "Teuchos_ParameterEntryXMLConverter.hpp"

namespace Teuchos{

ParameterEntry
ParameterEntryXMLConverter::fromXMLtoParameterEntry(const XMLObject &xmlObj) const
{
/*  TEST_FOR_EXCEPTION(xmlObj.getRequired(getTypeAttributeName()) != getTypeAttributeValue(), 
    std::runtime_error, 
    "This converter is not approriate for converting a ParameterEntry tag with a type of " 
    << xmlObj.getRequired(getTypeAttributeName()) << " to a ParameterEntry with type " 
    << getTypeAttributeValue());*/
  // 2010/07/30: rabartl: Remove dead code above.
  ParameterEntry toReturn;
  bool isDefault = false;
  bool isUsed = false;


  if(xmlObj.hasAttribute(getDefaultAttributeName())){
    isDefault = xmlObj.getRequiredBool(getDefaultAttributeName());
  }

  if(xmlObj.hasAttribute(getUsedAttributeName())){
    isUsed = xmlObj.getRequiredBool(getUsedAttributeName());
  }

  setEntryValue(toReturn, xmlObj, isDefault);
  
  if(isUsed){
    toReturn.getAny();
  }
  
  return toReturn;
}


XMLObject
ParameterEntryXMLConverter::fromParameterEntrytoXML(const ParameterEntry &entry,
  const std::string &name) const
{
  //TEST_FOR_EXCEPTION(!isAppropriateConverter(entry), std::runtime_error, "This converter is not approriate for converting the ParameterEntry " << name << " to the xml tag with a type attribute of " << getTypeAttributeValue());
  XMLObject toReturn(ParameterEntry::getTagName());
  toReturn.addAttribute(getNameAttributeName(), name);
  toReturn.addAttribute(getTypeAttributeName(), getTypeAttributeValue());
  toReturn.addAttribute(getValueAttributeName(), getValueAttributeValue(entry));
  toReturn.addBool(getDefaultAttributeName(), entry.isDefault());
  toReturn.addBool(getUsedAttributeName(), entry.isUsed());
  return toReturn;
}


} // namespace Teuchos

