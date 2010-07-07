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

#ifndef TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardParameterEntryXMLConverters.hpp
*/
namespace Teuchos {

template<class IntegralType>
class StringToIntegralValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
	std::string getTagName() const;
private:
	static std::string& tagName(){
		static std::string tagName;
		if(tagName == ""){
			tagName = "stringtointegralvalidator";
		}
		return tagName;
	}

	static std::string& stringsTagName(){
		static std::string tagName;
		if(tagName == ""){
			tagName = "strings";
		}
		return tagName;
	}

	static std::string& stringTagName(){
		static std::string tagName;
		if(tagName == ""){
			tagName = "string";
		}
		return tagName;
	}

	static std::string& stringDocAttributeName(){
		static std::string tagName;
		if(tagName == ""){
			tagName = "stringdoc";
		}
		return tagName;
	}

	static std::string& integralValueAttributeName(){
		static std::string tagName;
		if(tagName == ""){
			tagName = "integralvalue";
		}
		return tagName;
	}

	static std::string& defaultParameterAttributeName(){
		static std::string attributeName;
		if(attributeName == ""){
			attributeName = "defaultparametername";
		}
		return attributeName;
	}
};

class AnyNumberValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
	std::string getTagName() const;
private:
	static std::string& tagName(){
		static std::string tagName;
		if(tagName == ""){
			tagName = "anynumbervalidator"
		}
		return tagName;
};

template<class S>
class EnhancedNumberValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
	std::string getTagName() const;
};

}
#endif // TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
