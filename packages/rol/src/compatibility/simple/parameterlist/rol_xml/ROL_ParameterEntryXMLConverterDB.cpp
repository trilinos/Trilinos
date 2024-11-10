// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"

using namespace Teuchos;

namespace ROL {

/**
 * \brief Thrown when a parameter entry tag is missing it's type attribute.
 */
class NoTypeAttributeException : public std::logic_error{
public:
  /**
   * \brief Constructs a NoTypeAttributeException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoTypeAttributeException(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when an appropriate ParameterEntryXMLConverter
 * can't be found.
 */
class CantFindParameterEntryConverterException : public std::logic_error{

public:

  /**
   * \brief Constructs an CantFindParameterEntryConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindParameterEntryConverterException(const std::string& what_arg):std::logic_error(what_arg){}

};


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
  TEUCHOS_TEST_FOR_EXCEPTION(
    !xmlObject.hasAttribute(ParameterEntryXMLConverter::getTypeAttributeName()),
    NoTypeAttributeException,
    ParameterEntry::getTagName() <<" tags must "
    "have a " << ParameterEntryXMLConverter::getTypeAttributeName() <<
    " attribute." << std::endl <<
    "Bad Parameter: " <<
    xmlObject.getAttribute(ParameterEntryXMLConverter::getNameAttributeName()) <<
    std::endl << std::endl);

  std::string parameterType = xmlObject.getRequired(
    ParameterEntryXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(parameterType);

  TEUCHOS_TEST_FOR_EXCEPTION(it == getConverterMap().end(),
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


} //namespace ROL


namespace {

TEUCHOS_STATIC_SETUP()
{
    typedef unsigned int uint;
    typedef unsigned short int ushort;
    typedef unsigned long ulong;
    TEUCHOS_ADD_TYPE_CONVERTER(int);
    TEUCHOS_ADD_TYPE_CONVERTER(uint);
    TEUCHOS_ADD_TYPE_CONVERTER(short);
    TEUCHOS_ADD_TYPE_CONVERTER(ushort);
    TEUCHOS_ADD_TYPE_CONVERTER(long);
    TEUCHOS_ADD_TYPE_CONVERTER(ulong);
    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    typedef long long int llint;
    typedef unsigned long long int ullint;
    TEUCHOS_ADD_TYPE_CONVERTER(llint);
    TEUCHOS_ADD_TYPE_CONVERTER(ullint);
    #endif //HAVE_TEUCHOS_LONG_LONG_INT
    TEUCHOS_ADD_TYPE_CONVERTER(double);
    TEUCHOS_ADD_TYPE_CONVERTER(float);

    typedef std::string myString;
    TEUCHOS_ADD_TYPE_CONVERTER(myString);

    TEUCHOS_ADD_TYPE_CONVERTER(char);
    TEUCHOS_ADD_TYPE_CONVERTER(bool);

    ROL::ParameterEntryXMLConverterDB::addConverter(
      Teuchos::rcp(new ROL::AnyParameterEntryConverter));
}


} // namespace
