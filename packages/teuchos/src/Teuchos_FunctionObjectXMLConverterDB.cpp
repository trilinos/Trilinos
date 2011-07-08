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

#include "Teuchos_FunctionObjectXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"
#include "Teuchos_XMLFunctionObjectExceptions.hpp"



namespace Teuchos {


void FunctionObjectXMLConverterDB::addConverter(
  RCP<const FunctionObject> function,
  RCP<FunctionObjectXMLConverter> converterToAdd){
  getConverterMap().insert(
    ConverterPair(function->getTypeAttributeValue(), converterToAdd));
}


RCP<const FunctionObjectXMLConverter>
FunctionObjectXMLConverterDB::getConverter(const FunctionObject& function){
  ConverterMap::const_iterator it = 
    getConverterMap().find(function.getTypeAttributeValue());
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindFunctionObjectConverterException,
    "Could not find a FunctionObjectXMLConverter for a FuncitonObject of type " <<
    function.getTypeAttributeValue() << " when writing out a condition to " <<
    "xml." << std::endl << std::endl
  )
  return it->second;
}


RCP<const FunctionObjectXMLConverter>
FunctionObjectXMLConverterDB::getConverter(const XMLObject& xmlObject)
{ 
  std::string functionType = xmlObject.getRequired(
    FunctionObjectXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(functionType);
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindFunctionObjectConverterException,
    "Could not find a FunctionObjectXMLConverter for a condition of type " <<
    functionType << " when reading in a condition from " <<
    "xml." << std::endl << std::endl
  )
  return it->second;
}

XMLObject FunctionObjectXMLConverterDB::convertFunctionObject(
  RCP<const FunctionObject> function)
{
  return getConverter(*function)->fromFunctionObjecttoXML(function);
}
 
RCP<FunctionObject> FunctionObjectXMLConverterDB::convertXML(
  const XMLObject& xmlObject)
{
  return FunctionObjectXMLConverterDB::getConverter(xmlObject)->
    fromXMLtoFunctionObject(xmlObject);
}

FunctionObjectXMLConverterDB::ConverterMap&
FunctionObjectXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  return masterMap;
}


} // namespace Teuchos


namespace {


TEUCHOS_STATIC_SETUP()
{
    
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(short int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned short int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(long int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned long int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(double);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(float);

    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(long long int);
    TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(unsigned long long int);
    #endif // HAVE_TEUCHOS_LONG_LONG_INT

}


} // namespace
