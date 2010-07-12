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


#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardValidatorXMLConverters.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

/*! \file Teuchos_ValidatorXMLCoverterDB.hpp
*/

namespace Teuchos {
#define ADD_STRINGTOINTEGRALCONVERTER(INTEGRALTYPE, PREFIXNAME) \
	\
	StringToIntegralParameterEntryValidator< INTEGRALTYPE > sti##PREFIXNAME##Validator(dummyStringArray, dummyDefaultName); \
	masterMap.insert(ConverterPair(sti##PREFIXNAME##Validator.getXMLTagName(), rcp(new StringToIntegralValidatorXMLConverter< INTEGRALTYPE >)));

#define ADD_ENHANCEDNUMBERCONVERTER(T, PREFIXNAME) \
	\
	EnhancedNumberValidator< T > en##PREFIXNAME##Validator; \
	masterMap.insert(ConverterPair(en##PREFIXNAME##Validator.getXMLTagName(), rcp(new EnhancedNumberValidatorXMLConverter< T >)));

#define ADD_ARRAYCONVERTER(VALIDATORTYPE, ENTRYTYPE, PREFIXNAME) \
	\
	ArrayValidator< VALIDATORTYPE , ENTRYTYPE > array##PREFIXNAME##Validator(rcp(new VALIDATORTYPE)); \
	masterMap.insert(ConverterPair(array##PREFIXNAME##Validator.getXMLTagName(), rcp(new ArrayValidatorXMLConverter< VALIDATORTYPE, ENTRYTYPE >)));

	void ValidatorXMLConverterDB::addConverter(ParameterEntryValidator& validator, RCP<ValidatorXMLConverter> converterToAdd){
		getConverterMap().insert(ConverterPair(validator.getXMLTagName(), converterToAdd));
	}

	RCP<const ValidatorXMLConverter> ValidatorXMLConverterDB::getConverter(const ParameterEntryValidator& validator){
		ConverterMap::const_iterator it = getConverterMap().find(validator.getXMLTagName());
		if(it != getConverterMap().end()){
			return it->second;
		}
		return getDefaultConverter();
	}

	RCP<const ValidatorXMLConverter> ValidatorXMLConverterDB::getConverter(const XMLObject& xmlObject){ 
		std::string parameterType = xmlObject.getTag();
		ConverterMap::const_iterator it = getConverterMap().find(parameterType);
		if(it != getConverterMap().end()){
			return it->second;
		}
		return getDefaultConverter();
	}

	RCP<ValidatorXMLConverter> ValidatorXMLConverterDB::getDefaultConverter(){
		static RCP<ValidatorXMLConverter> defaultConverter;
		if(defaultConverter.is_null()){
			defaultConverter = rcp(new UnknownValidatorXMLConverter);
		}
		return defaultConverter;
	}

	ValidatorXMLConverterDB::ConverterMap& ValidatorXMLConverterDB::getConverterMap(){
		static ConverterMap masterMap;
		if(masterMap.size() == 0){
			std::string dummyDefaultName = "";
			Array<std::string> dummyStringArray;
			ADD_STRINGTOINTEGRALCONVERTER(int, Int);
			ADD_STRINGTOINTEGRALCONVERTER(unsigned int, UnsignedInt);
			ADD_STRINGTOINTEGRALCONVERTER(short int, Short);
			ADD_STRINGTOINTEGRALCONVERTER(unsigned short int, UnsignedShort);
			ADD_STRINGTOINTEGRALCONVERTER(long int, Long);
			ADD_STRINGTOINTEGRALCONVERTER(unsigned long int, UnsignedLong);
			ADD_STRINGTOINTEGRALCONVERTER(double, Double);
			ADD_STRINGTOINTEGRALCONVERTER(float, Float);

			ADD_ENHANCEDNUMBERCONVERTER(int, Int);
			ADD_ENHANCEDNUMBERCONVERTER(unsigned int, UnsignedInt);
			ADD_ENHANCEDNUMBERCONVERTER(short int, Short);
			ADD_ENHANCEDNUMBERCONVERTER(unsigned short int, UnsignedShort);
			ADD_ENHANCEDNUMBERCONVERTER(long int, Long);
			ADD_ENHANCEDNUMBERCONVERTER(unsigned long int, UnsignedLong);
			ADD_ENHANCEDNUMBERCONVERTER(double, Double);
			ADD_ENHANCEDNUMBERCONVERTER(float, Float);

			ADD_ARRAYCONVERTER(EnhancedNumberValidator<int>, int, Int);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<unsigned int>, unsigned int, UnsignedInt);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<short int>, short int, Short);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<unsigned short int>, unsigned short int, UnsignedShort);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<long int>, long int, Long);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<unsigned long int>, unsigned long int, UnsignedLong);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<double>, double, Double);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<float>, float, Float);

			ADD_ARRAYCONVERTER(FileNameValidator, std::string, fileName);
			ADD_ARRAYCONVERTER(StringValidator, std::string, string);

			#ifdef HAVE_TEUCHOS_LONG_LONG_INT
			ADD_STRINGTOINTEGRALCONVERTER(long long int, LongLong);
			ADD_STRINGTOINTEGRALCONVERTER(unsigned long long int, UnsignedLongLong);
			ADD_ENHANCEDNUMBERCONVERTER(long long int, LongLong);
			ADD_ENHANCEDNUMBERCONVERTER(unsigned long long int, UnsignedLongLong);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<long long int>, long long int, LongLong);
			ADD_ARRAYCONVERTER(EnhancedNumberValidator<unsigned long long int>, unsigned long long int, UnsignedLongLong);
			#endif // HAVE_TEUCHOS_LONG_LONG_INT
	
			FileNameValidator fileNameValidator;
			masterMap.insert(ConverterPair(fileNameValidator.getXMLTagName(), rcp(new FileNameValidatorXMLConverter)));

			StringValidator stringValidator;
			masterMap.insert(ConverterPair(stringValidator.getXMLTagName(), rcp(new StringValidatorXMLConverter)));

			AnyNumberParameterEntryValidator anyNumberValidator;
			masterMap.insert(ConverterPair(anyNumberValidator.getXMLTagName(), rcp(new AnyNumberValidatorXMLConverter)));
		}
		return masterMap;
	}

}// end namespace Teuchos

