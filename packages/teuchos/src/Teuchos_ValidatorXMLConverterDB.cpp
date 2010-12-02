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
#include "Teuchos_StaticSetupMacro.hpp"



namespace Teuchos {


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
  return masterMap;
  // See default setup code below!
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


namespace {


TEUCHOS_STATIC_SETUP()
{
  Teuchos::ValidatorXMLConverterDB::ConverterMap& masterMap = 
    Teuchos::ValidatorXMLConverterDB::getConverterMap();
  TEUCHOS_ADD_NUMBERTYPECONVERTERS(masterMap, int);
  TEUCHOS_ADD_ENHANCEDNUMBERCONVERTER(masterMap, double);
  TEUCHOS_ADD_ENHANCEDNUMBERCONVERTER(masterMap, float);
  
  TEUCHOS_ADD_ARRAYCONVERTER(masterMap, Teuchos::EnhancedNumberValidator<double>, double);
  TEUCHOS_ADD_ARRAYCONVERTER(masterMap, Teuchos::EnhancedNumberValidator<float>, float);
  
  TEUCHOS_ADD_ARRAYCONVERTER(masterMap, Teuchos::FileNameValidator, std::string);
  TEUCHOS_ADD_ARRAYCONVERTER(masterMap, Teuchos::StringValidator, std::string);
  
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  TEUCHOS_ADD_NUMBERTYPECONVERTERS(masterMap, long long int);
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  
  TEUCHOS_ADD_STRINGTOINTEGRALCONVERTER(masterMap, Teuchos::EVerbosityLevel); 

  TEUCHOS_ADD_VALIDATOR_CONVERTER(masterMap, Teuchos::FileNameValidator, Teuchos::FileNameValidatorXMLConverter);
  TEUCHOS_ADD_VALIDATOR_CONVERTER(masterMap, Teuchos::StringValidator, Teuchos::StringValidatorXMLConverter);
  TEUCHOS_ADD_VALIDATOR_CONVERTER(masterMap, Teuchos::AnyNumberParameterEntryValidator, Teuchos::AnyNumberValidatorXMLConverter);
  
}


} // namespace
