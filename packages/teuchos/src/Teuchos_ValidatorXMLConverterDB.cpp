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

#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_StaticSetupMacro.hpp"



namespace Teuchos {


void ValidatorXMLConverterDB::addConverter(
  RCP<const ParameterEntryValidator> validator,
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
  TEUCHOS_ADD_NUMBERTYPECONVERTERS(int);
  TEUCHOS_ADD_ENHANCEDNUMBERCONVERTER(double);
  TEUCHOS_ADD_ENHANCEDNUMBERCONVERTER(float);
  
  TEUCHOS_ADD_ARRAYCONVERTER(Teuchos::EnhancedNumberValidator<double>, double);
  TEUCHOS_ADD_ARRAYCONVERTER(Teuchos::EnhancedNumberValidator<float>, float);
  
  TEUCHOS_ADD_ARRAYCONVERTER(Teuchos::FileNameValidator, std::string);
  TEUCHOS_ADD_ARRAYCONVERTER(Teuchos::StringValidator, std::string);
  
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  TEUCHOS_ADD_NUMBERTYPECONVERTERS(long long int);
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  
  TEUCHOS_ADD_STRINGTOINTEGRALCONVERTER(Teuchos::EVerbosityLevel); 

  TEUCHOS_ADD_VALIDATOR_CONVERTER(Teuchos::FileNameValidator, Teuchos::FileNameValidatorXMLConverter);
  TEUCHOS_ADD_VALIDATOR_CONVERTER(Teuchos::StringValidator, Teuchos::StringValidatorXMLConverter);
  TEUCHOS_ADD_VALIDATOR_CONVERTER(Teuchos::AnyNumberParameterEntryValidator, Teuchos::AnyNumberValidatorXMLConverter);
  
}


} // namespace
