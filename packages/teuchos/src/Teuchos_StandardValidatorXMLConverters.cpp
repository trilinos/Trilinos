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
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos{

	RCP<ParameterEntryValidator> IntStringToIntegralValidatorXMLConverter::fromXMLtoValidator(const XMLObject& validator) const{
		return StringToIntegralValidatorXMLConverter::fromXMLtoValidator<int>(validator);
	}

	XMLObject IntStringToIntegralValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::fromValidatortoXML<int>(validator);
	}

	bool IntStringToIntegralValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::isAppropriateConverter<int>(validator);
	}

	RCP<ParameterEntryValidator> ShortStringToIntegralValidatorXMLConverter::fromXMLtoValidator(const XMLObject& validator) const{
		return StringToIntegralValidatorXMLConverter::fromXMLtoValidator<short>(validator);
	}

	XMLObject ShortStringToIntegralValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::fromValidatortoXML<short>(validator);
	}

	bool ShortStringToIntegralValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::isAppropriateConverter<short>(validator);
	}

	RCP<ParameterEntryValidator> LongStringToIntegralValidatorXMLConverter::fromXMLtoValidator(const XMLObject& validator) const{
		return StringToIntegralValidatorXMLConverter::fromXMLtoValidator<long>(validator);
	}

	XMLObject LongStringToIntegralValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::fromValidatortoXML<long>(validator);
	}

	bool LongStringToIntegralValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::isAppropriateConverter<long>(validator);
	}

	RCP<ParameterEntryValidator> DoubleStringToIntegralValidatorXMLConverter::fromXMLtoValidator(const XMLObject& validator) const{
		return StringToIntegralValidatorXMLConverter::fromXMLtoValidator<double>(validator);
	}

	XMLObject DoubleStringToIntegralValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::fromValidatortoXML<double>(validator);
	}

	bool DoubleStringToIntegralValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::isAppropriateConverter<double>(validator);
	}

	RCP<ParameterEntryValidator> FloatStringToIntegralValidatorXMLConverter::fromXMLtoValidator(const XMLObject& validator) const{
		return StringToIntegralValidatorXMLConverter::fromXMLtoValidator<float>(validator);
	}

	XMLObject FloatStringToIntegralValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::fromValidatortoXML<float>(validator);
	}

	bool FloatStringToIntegralValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return StringToIntegralValidatorXMLConverter::isAppropriateConverter<float>(validator);
	}


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

	XMLObject AnyNumberValidatorXMLConverter::fromValidatortoXML(const RCP<ParameterEntryValidator> validator) const{
		TEST_FOR_EXCEPTION(!isAppropriateConverter(validator), std::runtime_error, "An AnyNumberValidatorXMLConverter is not apporpriate for this type of validator.");
		RCP<AnyNumberParameterEntryValidator> convertedValidator = rcp_static_cast<AnyNumberParameterEntryValidator>(validator);
		XMLObject toReturn(validator->getXMLTagName());
		toReturn.addBool(getAllowIntAttributeName(), convertedValidator->allowInt());
		toReturn.addBool(getAllowDoubleAttributeName(), convertedValidator->allowDouble());
		toReturn.addBool(getAllowStringAttributeName(), convertedValidator->allowString());
		toReturn.addAttribute(getPrefferedTypeAttributeName(), convertedValidator->getPrefferedTypeString(convertedValidator->prefferedType()));
		return toReturn;
	}

	bool AnyNumberValidatorXMLConverter::isAppropriateConverter(const RCP<ParameterEntryValidator> validator) const{
		return !(rcp_dynamic_cast<AnyNumberParameterEntryValidator>(validator).is_null());
	}
}

