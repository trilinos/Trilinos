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
