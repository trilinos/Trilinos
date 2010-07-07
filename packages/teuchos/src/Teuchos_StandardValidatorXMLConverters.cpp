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
#include "Teuchos_StandardValidatorXMLConverters.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos{

	RCP<ParameterEntryValidator> StringToIntegralValidatorXMLConverter::fromXMLtoValidator(const XMLObject& xmlObj) const{
		TEST_FOR_EXCEPTION(xmlObj.getTag() != getTagName(), 
			std::runtime_error, 
			"Cannot convert xmlObject to StringToIntegralValidator. Expected a " << getTagName() 
			<< " tag but got a " << xmlObj.getTag() << "tag");
		TEST_FOR_EXCEPTION(xmlObj.child(0).getTag() != stringsTagName(), 
			std::runtime_error, 
			"Cannot convert xmlObject to StringToIntegralValidator. The " << getTagName() 
			<< " tag's first child should be a " << stringsTagName() << "tag");
		std::vector<std::string> strings;
		std::vector<std::string> stringDocs;

	}

	XMLObject StringToIntegralValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		RCP<StringToIntegralParameterEntryValidator<IntegralType> > convertedValidator =
			rcp_static_cast<StringToIntegralValidatorXMLConverter<IntegralType> >(validator);
		XMLObject toReturn(getTagName());
		XMLObject stringsTag(stringsTagName());
		Array<std::string>::const_iterator it = convertedValidator->validStringValues()->begin();
		for(;
			it != convertedValidator->validStringValues()->end(); ++it){
			XMLObject stringTag(stringTagName());
			stringTag.addContent(*it);
			stringTag.addAttribute(integralValueAttributeName(), convertedValidator->getIntegralValue(*it));
			stringsTag.addChild(stringTag);
		}
		if(convertedValidator->getStringDocs()->size()!=0){
			Array<std::string>::const_iterator it2 = convertedValidator->getStringDocs()->begin();
			for(int i=0; i<stringsTag.numChildren() && it2 != convertedValidator->getStringDocs()->end(); ++i; ++it2){
				stringsTag.getChild(i).addAttribute(stringDocAttributeName(), *it2);		
			}
		}
		toReturn.addAttribute(defaultParameterAttributeName(), convertedValidator->getDefaultParameterName());
		toReturn.addChild(stringsTag);
		return toReturn;
	}

	bool StringToIntegralValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<StringToIntegralValidator<IntegralType> >(validator).is_null()) 
	}

	std::string StringToIntegralValidatorXMLConverter::getTagName() const{
		return tagName();
	}

	RCP<ParameterEntryValidator> AnyNumberValidatorXMLConverter::fromXMLtoValidator(const XMLObject& xmlObj) const{
		TEST_FOR_EXCEPTION(xmlObj.getTag() != getTagName(), 
			std::runtime_error, 
			"Cannot convert xmlObject to StringToIntegralValidator. Expected a " << getTagName() 
			<< " tag but got a " << xmlObj.getTag() << "tag");

	}

	XMLObject AnyNumberValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		TEST_FOR_EXCEPTION(!isAppropriateConverter, std::runtime_error, "An AnyNumberValidatorXMLConverter is not apporpriate for this type of validator.");
		RCP<AnyNumberValidator> convertedValidator = rcp_static_cast<AnyNumberValidator>(validator);
		XMLObject toReturn(getTagName());
		return toReturn;
	}

	bool AnyNumberValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<AnyNumberValidator>(validator).is_null()) 
	}

	std::string AnyNumberValidatorXMLConverter::getTagName() const{
		return tagName();
	}
}

