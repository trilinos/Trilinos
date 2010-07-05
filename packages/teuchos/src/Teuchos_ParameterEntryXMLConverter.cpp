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

ParameterEntry ParameterEntryXMLConverter::fromXMLtoParameterEntry(const XMLObject &xmlObj) const{
	ParameterEntry toReturn;
	bool isDefault = false;
	bool isUsed = false;

	if(xmlObj.hasAttribute("isDefault")){
		isDefault = xmlObj.getRequiredBool("isDefault");
	}

	if(xmlObj.hasAttribute("isUsed")){
		isUsed = xmlObj.getRequiredBool("isUsed");
	}

	setEntryValue(toReturn, xmlObj, isDefault);
	
	if(isUsed){
		toReturn.getAny();
	}
	
	return toReturn;
}

XMLObject ParameterEntryXMLConverter::fromParameterEntrytoXML(const ParameterEntry &entry, const std::string &name) const{
	XMLObject toReturn(ParameterEntry::getTagName());
	toReturn.addAttribute("name", name);
	toReturn.addAttribute("type", getTypeAttributeValue());
	toReturn.addAttribute("value", getValueAttributeValue(entry));
	toReturn.addBool("isDefault",entry.isDefault());
	toReturn.addBool("isUsed", entry.isUsed());
	return toReturn;
}

}

