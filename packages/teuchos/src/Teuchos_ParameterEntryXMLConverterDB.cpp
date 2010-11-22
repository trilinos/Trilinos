
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

#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"

namespace Teuchos{


#define ADD_TYPE_CONVERTER(T) \
  \
  RCP<StandardTemplatedParameterConverter< T > > T##Converter = \
    rcp(new StandardTemplatedParameterConverter< T >); \
  masterMap.insert(ConverterPair(T##Converter->getTypeAttributeValue(), \
      T##Converter)); 

#define ADD_ARRAYTYPE_CONVERTER(T) \
  RCP<StandardTemplatedParameterConverter< Array< T > > > T##ArrayConverter = \
    rcp(new StandardTemplatedParameterConverter< Array< T > >); \
  masterMap.insert(ConverterPair(T##ArrayConverter->getTypeAttributeValue(), \
      T##ArrayConverter)); 

#define ADD_TYPE_AND_ARRAYTYPE_CONVERTER(T) \
  \
  ADD_TYPE_CONVERTER(T); \
  ADD_ARRAYTYPE_CONVERTER(T);

RCP<const ParameterEntryXMLConverter> 
ParameterEntryXMLConverterDB::getConverter(RCP<const ParameterEntry> entry) 
{
  ConverterMap::const_iterator it = 
  getConverterMap().find(entry->getAny().typeName());
  if(it == getConverterMap().end()){
    return getDefaultConverter();
  }
  else{
    return it->second;
  }
} 

RCP<const ParameterEntryXMLConverter> 
ParameterEntryXMLConverterDB::getConverter(const XMLObject& xmlObject)
{ 
  TEST_FOR_EXCEPTION(
    !xmlObject.hasAttribute(ParameterEntryXMLConverter::getTypeAttributeName()), 
    NoTypeAttributeExecption,
    ParameterEntry::getTagName() <<" tags must "
    "have a " << ParameterEntryXMLConverter::getTypeAttributeName() << 
    " attribute." << std::endl <<
    "Bad Parameter: " << 
    xmlObject.getAttribute(XMLParameterListWriter::getNameAttributeName()) <<
    std::endl << std::endl);
 
  std::string parameterType = xmlObject.getRequired(
    ParameterEntryXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(parameterType);

  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindParameterEntryConverterException,
    "Can't find converter for parameter entry of type: " <<
    xmlObject.getRequired(ParameterEntryXMLConverter::getTypeAttributeName()) <<
    std::endl << std::endl);

  return it->second;
}

void ParameterEntryXMLConverterDB::printKnownConverters(std::ostream& out){
  out << "Known ParameterEntryXMLConverters: " << std::endl;
  for(
    ConverterMap::const_iterator it = getConverterMap().begin();
    it != getConverterMap().end();
    ++it)
  {
    out << "\t" << it->first <<std::endl;
  }
}
  
RCP<const ParameterEntryXMLConverter> 
  ParameterEntryXMLConverterDB::getDefaultConverter()
{
  static RCP<const AnyParameterEntryConverter> defaultConverter;
  if(defaultConverter.is_null()){ 
      defaultConverter = rcp(new AnyParameterEntryConverter);
  }
  return defaultConverter;
}
ParameterEntryXMLConverterDB::ConverterMap& 
ParameterEntryXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  if(masterMap.size() == 0){
    typedef unsigned int uint;
    typedef unsigned short int ushort;
    typedef unsigned long ulong;
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(int);
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(uint);
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(short);
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(ushort);
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(long);
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(ulong);
    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    typedef long long int llint;
    typedef unsigned long long int ullint;
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(llint);
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(ullint);
    #endif //HAVE_TEUCHOS_LONG_LONG_INT
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(double);
    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(float);

    ADD_TYPE_AND_ARRAYTYPE_CONVERTER(string);

    ADD_TYPE_CONVERTER(char);
    ADD_TYPE_CONVERTER(bool);

    RCP<AnyParameterEntryConverter > anyConverter = 
      rcp(new AnyParameterEntryConverter); 
    masterMap.insert(ConverterPair(anyConverter->getTypeAttributeValue(), 
      anyConverter)); 

  }
  return masterMap;
}
  


} //namespace Teuchos
