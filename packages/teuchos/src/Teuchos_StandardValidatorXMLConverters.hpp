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


#ifndef TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDVALIDATORXMLCONVERTERSL_HPP

#include "Teuchos_ValidatorXMLConverter.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"

namespace Teuchos {

//class ParameterEntryValidator;

template<class IntegralType>
class StringToIntegralValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const{
		Array<std::string> strings;
		Array<std::string> stringDocs;
		Array<IntegralType> integralValues;
		for(int i=0; i<xmlObj.numChildren(); ++i){
			XMLObject currentChild = xmlObj.getChild(i);
			TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(), 
				std::runtime_error,  
				"Cannot convert xmlObject to StringToIntegralValidator." 
				<< "\n Unrecognized tag: " << currentChild.getTag());
			strings.append(currentChild.getRequired(getStringValueAttributeName()));
			integralValues.append(currentChild.getWithDefault<IntegralType>(getIntegralValueAttributeName(),(IntegralType)i));
			stringDocs.append(currentChild.getWithDefault<std::string>(getStringDocAttributeName(),""));
		}
		std::string defaultParameterName = xmlObj.getRequired(getDefaultParameterAttributeName());
		return rcp(new StringToIntegralParameterEntryValidator<IntegralType>(strings, stringDocs, integralValues, defaultParameterName));
	}

	XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		RCP<const StringToIntegralParameterEntryValidator<IntegralType> > convertedValidator =
			rcp_static_cast<const StringToIntegralParameterEntryValidator<IntegralType> >(validator);
		XMLObject toReturn(validator->getXMLTagName());
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
			toReturn.addChild(stringTag);
		}
		toReturn.addAttribute(getDefaultParameterAttributeName(), convertedValidator->getDefaultParameterName());
		toReturn.addAttribute(getIntegralValueAttributeName(), TypeNameTraits<IntegralType>::name());
		return toReturn;
	}

	bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const StringToIntegralParameterEntryValidator<IntegralType> >(validator).is_null());
	}

private:
	static const std::string& getIntegralValueAttributeName(){
		static const std::string integralValueAttributeName_ = "integralvalue";
		return integralValueAttributeName_;
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

class AnyNumberValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const;
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

template<class T>
class EnhancedNumberValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const{
		RCP<EnhancedNumberValidator<T> > toReturn = rcp(new EnhancedNumberValidator<T>());
		xmlObj.getWithDefault(getStepAttributeName(), EnhancedNumberTraits<T>::defaultStep()),
		xmlObj.getWithDefault(getPrecisionAttributeName(), EnhancedNumberTraits<T>::defaultPrecision());
		if(xmlObj.hasAttribute(getMinAttributeName())){
			toReturn->setMin(xmlObj.getRequired<T>(getMinAttributeName()));
		}
		if(xmlObj.hasAttribute(getMaxAttributeName())){
			toReturn->setMax(xmlObj.getRequired<T>(getMaxAttributeName()));
		}
		return toReturn;
	}

	XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		RCP<const EnhancedNumberValidator<T> > convertedValidator =
			rcp_static_cast<const EnhancedNumberValidator<T> >(validator);
		XMLObject toReturn(convertedValidator->getXMLTagName());
		if(convertedValidator->hasMin()){
			toReturn.addAttribute<T>(getMinAttributeName(), convertedValidator->getMin());
		}
		if(convertedValidator->hasMax()){
			toReturn.addAttribute<T>(getMaxAttributeName(), convertedValidator->getMin());
		}
		toReturn.addAttribute<T>(getStepAttributeName(), convertedValidator->getStep());
		toReturn.addAttribute<T>(getPrecisionAttributeName(), convertedValidator->getPrecision());
		return toReturn;
	}

	bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const EnhancedNumberValidator<T> >(validator).is_null());
	}
private:

	static const std::string& getMinAttributeName(){
		static const std::string minAttributeName = "min";
		return minAttributeName;
	}

	static const std::string& getMaxAttributeName(){
		static const std::string maxAttributeName = "max";
		return maxAttributeName;
	}

	static const std::string& getStepAttributeName(){
		static const std::string stepAttributeName = "step";
		return stepAttributeName;
	}

	static const std::string& getPrecisionAttributeName(){
		static const std::string precisionAttributeName = "precision";
		return precisionAttributeName;
	}
};

class FileNameValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const;
private:
	static const std::string& getFileMustExistAttributeName(){
		static const std::string fileMustExistAttributeName = "filemustexist";
		return fileMustExistAttributeName;
	}
};

class StringValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const;
private:
	static const std::string& getStringTagName(){
		static const std::string stringTagName = "string";
		return stringTagName;
	}
	static const std::string& getStringValueAttributeName(){
		static const std::string stringValueAttributeName = "stringvalue";
		return stringValueAttributeName;
	}
};

template<class ValidatorType, class EntryType>
class ArrayValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const{
		RCP<ValidatorType> prototypeValidator = rcp_static_cast<ValidatorType>(ValidatorXMLConverterDB::getConverter(xmlObj.getChild(0).getTag())->fromXMLtoValidator(xmlObj.getChild(0)));
		return rcp(new ArrayValidator<ValidatorType, EntryType>(prototypeValidator));
	}

	XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		XMLObject toReturn(validator->getXMLTagName());
		RCP<const ArrayValidator<ValidatorType, EntryType> > convertedValidator = 
			rcp_static_cast<const ArrayValidator<ValidatorType, EntryType> >(validator);
		toReturn.addChild(ValidatorXMLConverterDB::getConverter(*(convertedValidator->getPrototype()))->fromValidatortoXML(convertedValidator->getPrototype()));
		return toReturn;
	}

	bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const ArrayValidator<ValidatorType, EntryType> >(validator).is_null());
	}
};

class UnknownValidatorXMLConverter : public ValidatorXMLConverter{
public:
	RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const;
	XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const;
	bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const;
};


} // end namespace Teuchos


#endif  // TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP

