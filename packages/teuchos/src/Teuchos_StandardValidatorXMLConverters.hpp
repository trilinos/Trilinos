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

#include "Teuchos_ValidatorXMLConverter.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#ifndef TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardParameterEntryXMLConverters.hpp
*/
namespace Teuchos {

class StringToIntegralValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const=0;

	template<class IntegralType>
	static RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj){
		RCP<ParameterEntryValidator> dummyValidator = getDummyValidator<IntegralType>();
		TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator->getXMLTagName(), 
			std::runtime_error, 
			"Cannot convert xmlObject to StringToIntegralValidator. Expected a " << dummyValidator->getXMLTagName() 
			<< " tag but got a " << xmlObj.getTag() << "tag");
		TEST_FOR_EXCEPTION(xmlObj.getChild(0).getTag() != getStringsTagName(), 
			std::runtime_error,  
			"Cannot convert xmlObject to StringToIntegralValidator. The " << dummyValidator->getXMLTagName() 
			<< " tag's first child should be a " << getStringsTagName() << "tag");
		XMLObject stringsTag = xmlObj.getChild(0);
		Array<std::string> strings;
		Array<std::string> stringDocs;
		Array<IntegralType> integralValues;
		for(int i=0; i<stringsTag.numChildren(); ++i){
			strings.append(stringsTag.getChild(i).getRequired(getStringValueAttributeName()));
			integralValues.append(stringsTag.getChild(i).getRequired<IntegralType>(getIntegralValueAttributeName()));
			if(stringsTag.hasAttribute(getStringDocAttributeName())) {
				stringDocs.append(stringsTag.getChild(i).getRequired(getStringDocAttributeName()));
			}
		}
		std::string defaultParameterName = xmlObj.getRequired(getDefaultParameterAttributeName());
		return rcp(new StringToIntegralParameterEntryValidator<IntegralType>(strings, stringDocs, integralValues, defaultParameterName));
	}	

	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const=0;

	template<class IntegralType>
	static XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator){
		RCP<StringToIntegralParameterEntryValidator<IntegralType> > convertedValidator =
			rcp_static_cast<StringToIntegralParameterEntryValidator<IntegralType> >(validator);
		XMLObject toReturn(validator->getXMLTagName());
		XMLObject stringsTag(getStringsTagName());
		RCP<const Array<std::string> > stringValues = convertedValidator->validStringValues();
		RCP<const Array<std::string> > stringDocValues = convertedValidator->getStringDocs();
		bool hasStringDocs = stringDocValues->size() != 0;
		for(int i =0; i<stringValues->size(); ++i){
			XMLObject stringTag(getStringTagName());
			stringTag.addAttribute(getStringValueAttributeName(), (*stringValues)[i]);
			stringTag.addAttribute(getIntegralValueAttributeName(), convertedValidator->getIntegralValue((*stringValues)[i]));
			if(hasStringDocs){
				stringTag.addAttribute(getStringDocAttributeName(), (*stringValues)[i]);
			}
			stringsTag.addChild(stringTag);
		}
		toReturn.addAttribute(getDefaultParameterAttributeName(), convertedValidator->getDefaultParameterName());
		toReturn.addAttribute(getIntegralValueAttributeName(), TypeNameTraits<IntegralType>::name());
		toReturn.addChild(stringsTag);
		return toReturn;
	}

	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const=0;

	template<class IntegralType>
	static bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator){
		return !(rcp_dynamic_cast<StringToIntegralParameterEntryValidator<IntegralType> >(validator).is_null());
	}	

	template<class IntegralType>
	static RCP<ParameterEntryValidator> getDummyValidator(){
		ArrayView<std::string> dummyStrings;
		std::string dummyParamName;
		return rcp(new StringToIntegralParameterEntryValidator<IntegralType>(dummyStrings, dummyParamName));
	}

private:
	static const std::string& getIntegralValueAttributeName(){
		static const std::string integralValueAttributeName_ = "integralvalue";
		return integralValueAttributeName_;
	}

	static const std::string& getStringsTagName(){
		static const std::string stringsTagName_ = "strings";
		return stringsTagName_;
	}

	static const std::string& getStringTagName(){
		static const std::string stringTagName_ = "string";
		return stringTagName_;
	}

	static const std::string& getStringValueAttributeName(){
		static const std::string stringValueAttributeName_ = "stringvalue";
		return stringValueAttributeName_;
	}

	static const std::string& getStringDocAttributeName(){
		static const std::string stringDocAttributeName_ = "stringdoc";
		return stringDocAttributeName_;
	}

	static const std::string& getDefaultParameterAttributeName(){
		static const std::string defaultParameterAttributeName_ = "defaultparametername";
		return defaultParameterAttributeName_;
	}
};

class IntStringToIntegralValidatorXMLConverter : public StringToIntegralValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
};

class ShortStringToIntegralValidatorXMLConverter : public StringToIntegralValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
};

class LongStringToIntegralValidatorXMLConverter : public StringToIntegralValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
};

class DoubleStringToIntegralValidatorXMLConverter : public StringToIntegralValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
};

class FloatStringToIntegralValidatorXMLConverter : public StringToIntegralValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
};

class AnyNumberValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
private:
	static const std::string& getAllowIntAttributeName(){
		static const std::string allowIntAttributeName_ = "allowInt";
		return allowIntAttributeName_;
	}

	static const std::string& getAllowDoubleAttributeName(){
		static const std::string allowDoubleAttributeName_ = "allowDouble";
		return allowDoubleAttributeName_;
	}

	static const std::string& getAllowStringAttributeName(){
		static const std::string allowStringAttributeName_ = "allowString";
		return allowStringAttributeName_;
	}
	
	static const std::string& getPrefferedTypeAttributeName(){
		static const std::string prefferedTypeAttributeName_ = "prefferedType";
		return prefferedTypeAttributeName_;
	}
};

template<class S>
class EnhancedNumberValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const;
};

}
#endif // TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
