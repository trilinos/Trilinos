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
#include "Teuchos_VerbosityLevel.hpp"



namespace Teuchos {

#define ADD_NUMBERTYPECONVERTERS(T) \
  ADD_STRINGTOINTEGRALCONVERTER( T ); \
  ADD_ENHANCEDNUMBERCONVERTER( T ); \
  ADD_ARRAYCONVERTER(EnhancedNumberValidator< T >, T );

#define ADD_STRINGTOINTEGRALCONVERTER(INTEGRALTYPE) \
  \
  masterMap.insert(ConverterPair( \
    DummyObjectGetter<StringToIntegralParameterEntryValidator< INTEGRALTYPE > >:: \
      getDummyObject()->getXMLTypeName(), \
    rcp(new StringToIntegralValidatorXMLConverter< INTEGRALTYPE >)));


#define ADD_ENHANCEDNUMBERCONVERTER(T) \
  \
  masterMap.insert(ConverterPair( \
    DummyObjectGetter<EnhancedNumberValidator< T > >:: \
      getDummyObject()->getXMLTypeName(), \
    rcp(new EnhancedNumberValidatorXMLConverter< T >)));


#define ADD_ARRAYCONVERTER( VALIDATORTYPE, ENTRYTYPE ) \
  \
  masterMap.insert(ConverterPair( \
    DummyObjectGetter<ArrayValidator< VALIDATORTYPE , ENTRYTYPE > >:: \
      getDummyObject()->getXMLTypeName(), \
    rcp(new ArrayValidatorXMLConverter< VALIDATORTYPE, ENTRYTYPE >)));

void ValidatorXMLConverterDB::addConverter(
  RCP<ParameterEntryValidator> validator,
  RCP<ValidatorXMLConverter> converterToAdd){
  getConverterMap().insert(ConverterPair(
    validator->getXMLTypeName(), converterToAdd));
}


RCP<const ValidatorXMLConverter>
ValidatorXMLConverterDB::getConverter(const ParameterEntryValidator& validator)
{
  ConverterMap::const_iterator it = getConverterMap().find(validator.getXMLTypeName());
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindValidatorConverterException,
    "Could not find a ValidatorXMLConverter for validator type " <<
     validator.getXMLTypeName() << std::endl <<
     "Try adding an appropriate converter to the ValidatorXMLConverterDB " <<
     "in order solve this problem." << std::endl << std::endl
  )
  return it->second;
}


RCP<const ValidatorXMLConverter>
ValidatorXMLConverterDB::getConverter(const XMLObject& xmlObject)
{ 
  std::string validatorType = xmlObject.getRequired(
    ValidatorXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(validatorType);
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindValidatorConverterException,
    "Could not find a ValidatorXMLConverter for type " << validatorType <<
    std::endl << 
    "Try adding an appropriate converter to the ValidatorXMLConverterDB " <<
    "in order solve this problem." << std::endl << std::endl
  )
  return it->second;
}

XMLObject ValidatorXMLConverterDB::convertValidator(
  RCP<const ParameterEntryValidator> validator,
  const ValidatortoIDMap& validatorIDsMap,
  bool assignID)
{
  return getConverter(*validator)->fromValidatortoXML(
    validator, validatorIDsMap, assignID);
}
 
RCP<ParameterEntryValidator> ValidatorXMLConverterDB::convertXML(
  const XMLObject& xmlObject,
  const IDtoValidatorMap& validatorIDsMap)
{
  return ValidatorXMLConverterDB::
    getConverter(xmlObject)->fromXMLtoValidator(xmlObject, validatorIDsMap);
}

ValidatorXMLConverterDB::ConverterMap&
ValidatorXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  if(masterMap.size() == 0){
    ADD_NUMBERTYPECONVERTERS(int);
    ADD_ENHANCEDNUMBERCONVERTER(double);
    ADD_ENHANCEDNUMBERCONVERTER(float);

    ADD_ARRAYCONVERTER(EnhancedNumberValidator<double>, double);
    ADD_ARRAYCONVERTER(EnhancedNumberValidator<float>, float);

    ADD_ARRAYCONVERTER(FileNameValidator, std::string);
    ADD_ARRAYCONVERTER(StringValidator, std::string);

    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    ADD_NUMBERTYPECONVERTERS(long long int);
    #endif // HAVE_TEUCHOS_LONG_LONG_INT

    ADD_STRINGTOINTEGRALCONVERTER( EVerbosityLevel ); \


    masterMap.insert(
      ConverterPair(
        DummyObjectGetter<FileNameValidator>::getDummyObject()->getXMLTypeName(), 
        rcp(new FileNameValidatorXMLConverter)));

    masterMap.insert(
      ConverterPair(
        DummyObjectGetter<StringValidator>::getDummyObject()->getXMLTypeName(), 
        rcp(new StringValidatorXMLConverter)));

    masterMap.insert(
      ConverterPair(
        DummyObjectGetter<AnyNumberParameterEntryValidator>::getDummyObject()->getXMLTypeName(), 
        rcp(new AnyNumberValidatorXMLConverter)));

  }
  return masterMap;
}

void ValidatorXMLConverterDB::printKnownConverters(std::ostream& out){
  out << "Known ValidatorXMLConverters: " << std::endl;
  for(
    ConverterMap::const_iterator it = getConverterMap().begin();
    it != getConverterMap().end();
    ++it)
  {
    out << "\t" << it->first <<std::endl;
  }
}
  

} // namespace Teuchos
