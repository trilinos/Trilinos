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

#include "Teuchos_DependencyXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"



namespace Teuchos {

  

void DependencyXMLConverterDB::addConverter(
  RCP<const Dependency> dependency,
  RCP<DependencyXMLConverter> converterToAdd)
{
  getConverterMap().insert(
    ConverterPair(dependency->getTypeAttributeValue(), converterToAdd));
}


RCP<const DependencyXMLConverter>
DependencyXMLConverterDB::getConverter(const Dependency& dependency)
{
  ConverterMap::const_iterator it = 
    getConverterMap().find(dependency.getTypeAttributeValue());
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindDependencyConverterException,
    "Could not find a DependencyXMLConverter for a dependency with "
    "attribute tag " << dependency.getTypeAttributeValue() << 
    "!" << std::endl <<
    "Try adding an appropriate converter to the DependencyXMLConverterDB " <<
    "in order to solve this problem." << std::endl << std::endl);
  return it->second;
}


RCP<const DependencyXMLConverter>
DependencyXMLConverterDB::getConverter(const XMLObject& xmlObject)
{ 
  std::string dependencyType = xmlObject.getRequired(
    DependencyXMLConverter::getTypeAttributeName());
  ConverterMap::const_iterator it = getConverterMap().find(dependencyType);
  TEST_FOR_EXCEPTION(it == getConverterMap().end(),
    CantFindDependencyConverterException,
    "Could not find a DependencyXMLConverter for a dependency of type " <<
    dependencyType << "!" << std::endl <<
    "Try adding an appropriate converter to the DependencyXMLConverterDB " <<
    "in order to solve this problem." << std::endl << std::endl);
  return it->second;
}

XMLObject DependencyXMLConverterDB::convertDependency(
  RCP<const Dependency> dependency,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
  ValidatortoIDMap& validatorIDsMap)
{
  return getConverter(*dependency)->fromDependencytoXML(
    dependency, entryIDsMap, validatorIDsMap);
}
 
RCP<Dependency> DependencyXMLConverterDB::convertXML(
    const XMLObject& xmlObject, 
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap)
{
  return DependencyXMLConverterDB::getConverter(xmlObject)->
    fromXMLtoDependency(xmlObject, entryIDsMap, validatorIDsMap);
}

DependencyXMLConverterDB::ConverterMap&
DependencyXMLConverterDB::getConverterMap()
{
  static ConverterMap masterMap;
  return masterMap;
}


} // namespace Teuchos


namespace {
  

TEUCHOS_STATIC_SETUP()
{
    TEUCHOS_ADD_TEMPLATED_NUMBER_DEPS(int)
    TEUCHOS_ADD_NUMBER_VISUAL_DEP(float)
    TEUCHOS_ADD_RANGE_VALIDATOR_DEP(float)
    TEUCHOS_ADD_NUMBER_VISUAL_DEP(double)
    TEUCHOS_ADD_RANGE_VALIDATOR_DEP(double)

    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    TEUCHOS_ADD_TEMPLATED_NUMBER_DEPS(long long int)
    #endif // HAVE_TEUCHOS_LONG_LONG_INT

    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::StringValidatorDependency,
      Teuchos::StringValidatorDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::StringVisualDependency,
      Teuchos::StringVisualDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::BoolValidatorDependency,
      Teuchos::BoolValidatorDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::BoolVisualDependency,
      Teuchos::BoolVisualDependencyXMLConverter)
    TEUCHOS_ADD_DEP_CONVERTER(
      Teuchos::ConditionVisualDependency,
      Teuchos::ConditionVisualDependencyXMLConverter)
}


} //namespace 

