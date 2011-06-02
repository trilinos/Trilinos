
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_StaticSetupMacro.hpp"
#include "Teuchos_TwoDArray.hpp"

namespace Teuchos{


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
  return masterMap;

}
  

} //namespace Teuchos


namespace {

TEUCHOS_STATIC_SETUP()
{
    typedef unsigned int uint;
    typedef unsigned short int ushort;
    typedef unsigned long ulong;
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(int);
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(uint);
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(short);
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(ushort);
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(long);
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(ulong);
    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    typedef long long int llint;
    typedef unsigned long long int ullint;
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(llint);
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(ullint);
    #endif //HAVE_TEUCHOS_LONG_LONG_INT
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(double);
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(float);

    typedef std::string myString;
    TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(myString);

    TEUCHOS_ADD_TYPE_CONVERTER(char);
    TEUCHOS_ADD_TYPE_CONVERTER(bool);

    Teuchos::ParameterEntryXMLConverterDB::addConverter(
      Teuchos::rcp(new Teuchos::AnyParameterEntryConverter));
}


} // namespace
