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

/*! \file Teuchos_ParameterEntryXMLCoverterDB.hpp
*/

namespace Teuchos {

#define ADD_TYPE_CONVERTER(T,PREFIXNAME) \
	\
	RCP<StandardTemplatedParameterConverter< T > > PREFIXNAME##Converter = rcp(new StandardTemplatedParameterConverter< T >); \
	masterMap.insert(ConverterPair(PREFIXNAME##Converter->getTypeAttributeValue(), PREFIXNAME##Converter)); 

#define ADD_ARRAYTYPE_CONVERTER(T,PREFIXNAME) \
	RCP<ArrayTemplatedParameterConverter< T > > PREFIXNAME##ArrayConverter = rcp(new ArrayTemplatedParameterConverter< T >); \
	masterMap.insert(ConverterPair(PREFIXNAME##ArrayConverter->getTypeAttributeValue(), PREFIXNAME##ArrayConverter)); 

#define ADD_TYPE_AND_ARRAYTYPE_CONVERTER(T , PREFIXNAME) \
	\
	ADD_TYPE_CONVERTER(T, PREFIXNAME); \
	ADD_ARRAYTYPE_CONVERTER(T,PREFIXNAME);

/* \class Teuchos::ParameterEntryXMLConverterDB
 *
 * \brief Provides ability to lookup ParameterEntryXMLConverters
 */
class ParameterEntryXMLConverterDB{
public:

	/* \brief Add a converter to the database.
	 *
	 * @param convertToAdd The converter to add to the database.
	 */
	static void addConverter(RCP<ParameterEntryXMLConverter> converterToAdd){
		getConverterMap().insert(ConverterPair(converterToAdd->getTypeAttributeValue(), converterToAdd));
	}

	/* \brief Get an appropriate ParameterEntryXMLConverter given a ParameterEntry.
	 *
	 * @param entry The ParameterEntry for which a converter is desired.
	 */
	static RCP<const ParameterEntryXMLConverter> getConverter(const ParameterEntry& entry){
		ConverterMap::const_iterator it = getConverterMap().find(entry.getAny().typeName());
		if(it != getConverterMap().end()){
			return it->second;
		}
		return getDefaultConverter();
	}

	/* \brief Get an appropriate ParameterEntryXMLConverter given a XMLObject.
	 *
	 * @param xmlObject The XMLObject for which a converter is desired.
	 */
	static RCP<const ParameterEntryXMLConverter> getConverter(const XMLObject& xmlObject){ 
		std::string parameterType = xmlObject.getRequired(ParameterEntryXMLConverter::getTypeAttributeName());
		ConverterMap::const_iterator it = getConverterMap().find(parameterType);
		if(it != getConverterMap().end()){
			return it->second;
		}
		return getDefaultConverter();
	}

private:
	/* \brief convience typedef */
	typedef std::map<std::string, RCP<ParameterEntryXMLConverter> > ConverterMap;
	/* \brief convience typedef */
	typedef std::pair<std::string, RCP<ParameterEntryXMLConverter> > ConverterPair;

	/* \brief Gets the default converter to be used to convert ParameterEntries. */
	static RCP<ParameterEntryXMLConverter> getDefaultConverter(){
		static RCP<ParameterEntryXMLConverter> defaultConverter;
		if(defaultConverter.is_null()){
			defaultConverter = rcp(new AnyParameterEntryConverter());
		}
		return defaultConverter;
	}

	/* \brief Gets the map containing all the converters. */
	static ConverterMap& getConverterMap(){
		static ConverterMap masterMap;
		if(masterMap.size() == 0){

			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(int, int);
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(unsigned int, unsignedInt);
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(short int, short);
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(unsigned short int, unsignedShortInt);
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(long int, long);
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(unsigned long int, unsignedLongInt);
			#ifdef HAVE_TEUCHOS_LONG_LONG_INT
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(long long int, longlong);
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(unsigned long long int, unsignedShortInt);
			#endif //HAVE_TEUCHOS_LONG_LONG_INT
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(double, double);
			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(float, float);

			ADD_TYPE_AND_ARRAYTYPE_CONVERTER(string, string);

			ADD_TYPE_CONVERTER(char, char);
			ADD_TYPE_CONVERTER(bool, bool);

		}
		return masterMap;
	}
};

}

#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
