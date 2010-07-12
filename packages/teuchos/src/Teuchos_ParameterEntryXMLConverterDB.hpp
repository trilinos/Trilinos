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


#ifndef TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
#define TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
#include "Teuchos_StandardParameterEntryXMLConverters.hpp"

/*! \file Teuchos_ParameterEntryXMLCoverter.hpp
*/

namespace Teuchos {

class ParameterEntryXMLConverterDB{
public:
	static void addConverter(RCP<ParameterEntryXMLConverter> converterToAdd){
		getConverterMap().insert(ConverterPair(converterToAdd->getTypeAttributeValue(), converterToAdd));
	}

	static RCP<const ParameterEntryXMLConverter> getConverter(const ParameterEntry& entry){
		ConverterMap::const_iterator it = getConverterMap().find(entry.getAny().typeName());
		if(it != getConverterMap().end()){
			return it->second;
		}
		return getDefaultConverter();
	}

	static RCP<const ParameterEntryXMLConverter> getConverter(const XMLObject& xmlObject){ 
		std::string parameterType = xmlObject.getRequired(ParameterEntryXMLConverter::getTypeAttributeName());
		ConverterMap::const_iterator it = getConverterMap().find(parameterType);
		if(it != getConverterMap().end()){
			return it->second;
		}
		return getDefaultConverter();
	}

private:
	typedef std::map<std::string, RCP<ParameterEntryXMLConverter> > ConverterMap;
	typedef std::pair<std::string, RCP<ParameterEntryXMLConverter> > ConverterPair;

	static RCP<ParameterEntryXMLConverter> getDefaultConverter(){
		static RCP<ParameterEntryXMLConverter> defaultConverter;
		if(defaultConverter.is_null()){
			defaultConverter = rcp(new AnyParameterEntryConverter());
		}
		return defaultConverter;
	}

	static ConverterMap& getConverterMap(){
		static ConverterMap masterMap;
		if(masterMap.size() == 0){
			RCP<IntParameterEntryConverter> intConverter = rcp(new IntParameterEntryConverter());
			masterMap.insert(ConverterPair(intConverter->getTypeAttributeValue(), intConverter));
			RCP<ShortParameterEntryConverter> shortConverter = rcp(new ShortParameterEntryConverter());
			masterMap.insert(ConverterPair(shortConverter->getTypeAttributeValue(), shortConverter));
			RCP<DoubleParameterEntryConverter> doubleConverter = rcp(new DoubleParameterEntryConverter());
			masterMap.insert(ConverterPair(doubleConverter->getTypeAttributeValue(), doubleConverter));
			RCP<FloatParameterEntryConverter> floatConverter = rcp(new FloatParameterEntryConverter());
			masterMap.insert(ConverterPair(floatConverter->getTypeAttributeValue(), floatConverter));
			RCP<StringParameterEntryConverter> stringConverter = rcp(new StringParameterEntryConverter());
			masterMap.insert(ConverterPair(stringConverter->getTypeAttributeValue(), stringConverter));
			RCP<CharParameterEntryConverter> charConverter = rcp(new CharParameterEntryConverter());
			masterMap.insert(ConverterPair(charConverter->getTypeAttributeValue(), charConverter));
			RCP<BoolParameterEntryConverter> boolConverter = rcp(new BoolParameterEntryConverter());
			masterMap.insert(ConverterPair(boolConverter->getTypeAttributeValue(), boolConverter));
			RCP<ArrayIntParameterEntryConverter> arrayintConverter = rcp(new ArrayIntParameterEntryConverter());
			masterMap.insert(ConverterPair(arrayintConverter->getTypeAttributeValue(), arrayintConverter));
			RCP<ArrayShortParameterEntryConverter> arrayshortConverter = rcp(new ArrayShortParameterEntryConverter());
			masterMap.insert(ConverterPair(arrayshortConverter->getTypeAttributeValue(), arrayshortConverter));
			RCP<ArrayDoubleParameterEntryConverter> arraydoubleConverter = rcp(new ArrayDoubleParameterEntryConverter());
			masterMap.insert(ConverterPair(arraydoubleConverter->getTypeAttributeValue(), arraydoubleConverter));
			RCP<ArrayFloatParameterEntryConverter> arrayfloatConverter = rcp(new ArrayFloatParameterEntryConverter());
			masterMap.insert(ConverterPair(arrayfloatConverter->getTypeAttributeValue(), arrayfloatConverter));
			RCP<ArrayStringParameterEntryConverter> arraystringConverter = rcp(new ArrayStringParameterEntryConverter());
			masterMap.insert(ConverterPair(arraystringConverter->getTypeAttributeValue(), arraystringConverter));
		}
		return masterMap;
	}
};

}

#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
