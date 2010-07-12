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
//#include "Teuchos_StandardParameterEntryValidators.hpp"
//#include "Teuchos_ValidatorXMLConverterDB.hpp"


/*! \file Teuchos_StandardValidatorXMLConverters.hpp
*/
namespace Teuchos {

/*	template<class IntegralType>
	RCP<ParameterEntryValidator> StringToIntegralValidatorXMLConverter<IntegralType>::fromXMLtoValidator(const XMLObject& xmlObj) const{
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

	template<class IntegralType>
	XMLObject StringToIntegralValidatorXMLConverter<IntegralType>::fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
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

	template<class IntegralType>
	bool StringToIntegralValidatorXMLConverter<IntegralType>::isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const StringToIntegralParameterEntryValidator<IntegralType> >(validator).is_null());
	}	
*/
	RCP<ParameterEntryValidator> AnyNumberValidatorXMLConverter::fromXMLtoValidator(const XMLObject& xmlObj) const{
		AnyNumberParameterEntryValidator dummyValidator;
		TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator.getXMLTagName(), 
			std::runtime_error, 
			"Cannot convert xmlObject to StringToIntegralValidator. Expected a " << dummyValidator.getXMLTagName() 
			<< " tag but got a " << xmlObj.getTag() << "tag");
		AnyNumberParameterEntryValidator::AcceptedTypes acceptedTypes;
		acceptedTypes.allowInt(xmlObj.getRequiredBool(getAllowIntAttributeName()));
		acceptedTypes.allowDouble(xmlObj.getRequiredBool(getAllowDoubleAttributeName()));
		acceptedTypes.allowString(xmlObj.getRequiredBool(getAllowStringAttributeName()));
		return rcp(new AnyNumberParameterEntryValidator(AnyNumberParameterEntryValidator::getPrefferedTypeStringEnum(xmlObj.getRequired(getPrefferedTypeAttributeName())), acceptedTypes));
	}

	XMLObject AnyNumberValidatorXMLConverter::fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		TEST_FOR_EXCEPTION(!isAppropriateConverter(validator), std::runtime_error, "An AnyNumberValidatorXMLConverter is not apporpriate for this type of validator.");
		RCP<const AnyNumberParameterEntryValidator> convertedValidator = rcp_static_cast<const AnyNumberParameterEntryValidator>(validator);
		XMLObject toReturn(validator->getXMLTagName());
		toReturn.addBool(getAllowIntAttributeName(), convertedValidator->allowInt());
		toReturn.addBool(getAllowDoubleAttributeName(), convertedValidator->allowDouble());
		toReturn.addBool(getAllowStringAttributeName(), convertedValidator->allowString());
		toReturn.addAttribute(getPrefferedTypeAttributeName(), convertedValidator->getPrefferedTypeString(convertedValidator->prefferedType()));
		return toReturn;
	}

	bool AnyNumberValidatorXMLConverter::isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(validator).is_null());
	}
/*
	template<class T>
	RCP<ParameterEntryValidator> EnhancedNumberValidatorXMLConverter<T>::fromXMLtoValidator(const XMLObject& xmlObj) const{
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

	template<class T>
	XMLObject EnhancedNumberValidatorXMLConverter<T>::fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
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

	template<class T>
	bool EnhancedNumberValidatorXMLConverter<T>::isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const EnhancedNumberValidator<T> >(validator).is_null());
	}	
*/
	RCP<ParameterEntryValidator> FileNameValidatorXMLConverter::fromXMLtoValidator(const XMLObject& xmlObj) const{
		FileNameValidator dummyValidator;
		TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator.getXMLTagName(), 
			std::runtime_error, 
			"Cannot convert xmlObject to StringToIntegralValidator. Expected a " << dummyValidator.getXMLTagName() 
			<< " tag but got a " << xmlObj.getTag() << "tag");
		return rcp(new FileNameValidator(xmlObj.getWithDefault<bool>(getFileMustExistAttributeName(), FileNameValidator::mustAlreadyExistDefault)));
	}

	XMLObject FileNameValidatorXMLConverter::fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		TEST_FOR_EXCEPTION(!isAppropriateConverter(validator), std::runtime_error, "An FileNameValidatorXMLConverter is not apporpriate for this type of validator.");
		RCP<const FileNameValidator> convertedValidator = rcp_static_cast<const FileNameValidator>(validator);
		XMLObject toReturn(validator->getXMLTagName());
		toReturn.addBool(getFileMustExistAttributeName(), convertedValidator->fileMustExist());
		return toReturn;
	}

	bool FileNameValidatorXMLConverter::isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const FileNameValidator>(validator).is_null());
	}

	RCP<ParameterEntryValidator> StringValidatorXMLConverter::fromXMLtoValidator(const XMLObject& xmlObj) const{
		StringValidator dummyValidator;
		TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator.getXMLTagName(), 
			std::runtime_error, 
			"Cannot convert xmlObject to StringValidator. Expected a " << dummyValidator.getXMLTagName() 
			<< " tag but got a " << xmlObj.getTag() << "tag");
		if(xmlObj.numChildren()!=0){
			Array<std::string> strings(xmlObj.numChildren());
			for(int i=0; i<xmlObj.numChildren(); ++i){
				XMLObject currentChild = xmlObj.getChild(i);
				TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(), 
					std::runtime_error,  
					"Cannot convert xmlObject to StringToIntegralValidator." 
					<< "\n Unrecognized tag: " << currentChild.getTag());
				strings[i] = (currentChild.getRequired(getStringValueAttributeName()));
			}
			return rcp(new StringValidator(strings));
		}
		return rcp(new StringValidator());
	}

	XMLObject StringValidatorXMLConverter::fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		TEST_FOR_EXCEPTION(!isAppropriateConverter(validator), std::runtime_error, "An StringValidatorXMLConverter is not apporpriate for this type of validator.");
		XMLObject toReturn(validator->getXMLTagName());
		Array<std::string>::const_iterator it = validator->validStringValues()->begin();
		for(; it != validator->validStringValues()->end(); ++it){
			XMLObject stringTag(getStringTagName());
			stringTag.addAttribute(getStringValueAttributeName(), *it);
			toReturn.addChild(stringTag);
		}
		return toReturn;
	}

	bool StringValidatorXMLConverter::isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const StringValidator>(validator).is_null());
	}

	RCP<ParameterEntryValidator> UnknownValidatorXMLConverter::fromXMLtoValidator(const XMLObject& xmlObj) const{
		throw std::runtime_error("Unknown xml tag. Can't convert to a Validator.");
		return Teuchos::null;
	}

	XMLObject UnknownValidatorXMLConverter::fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		throw std::runtime_error("Unknown validator. Convert to XML.");
		return NULL;
	}

	bool UnknownValidatorXMLConverter::isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return true;
	}
/*
	template<class ValidatorType, class EntryType>
	RCP<ParameterEntryValidator> ArrayValidatorXMLConverter<ValidatorType, EntryType>::fromXMLtoValidator(const XMLObject& xmlObj) const{
		RCP<ValidatorType> prototypeValidator = rcp_static_cast<ValidatorType>(ValidatorXMLConverterDB::getConverter(xmlObj.getChild(0).getTag())->fromXMLtoValidator(xmlObj.getChild(0)));
		return rcp(new ArrayValidator<ValidatorType, EntryType>(prototypeValidator));
	}

	template<class ValidatorType, class EntryType>
	XMLObject ArrayValidatorXMLConverter<ValidatorType, EntryType>::fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const{
		XMLObject toReturn(validator->getXMLTagName());
		RCP<const ArrayValidator<ValidatorType, EntryType> > convertedValidator = 
			rcp_static_cast<const ArrayValidator<ValidatorType, EntryType> >(validator);
		toReturn.addChild(ValidatorXMLConverterDB::getConverter(*(convertedValidator->getPrototype()))->fromValidatortoXML(convertedValidator->getPrototype()));
		return toReturn;
	}

	template<class ValidatorType, class EntryType>
	bool ArrayValidatorXMLConverter<ValidatorType, EntryType>::isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<const ArrayValidator<ValidatorType, EntryType> >(validator).is_null());
	}	*/

}

