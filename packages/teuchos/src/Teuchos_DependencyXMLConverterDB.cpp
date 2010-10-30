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

#include "Teuchos_DependencyXMLConverterDB.hpp"
#include "Teuchos_StandardDependencyXMLConverters.hpp"
#include "Teuchos_StandardDependencies.hpp"



namespace Teuchos {

#define ADD_DEP_TO_MAP(DEP_TYPE) \
    masterMap.insert( \
      ConverterPair( \
        DummyObjectGetter<DEP_TYPE>:: \
          getDummyObject()->getTypeAttributeValue(), \
        rcp(new DEP_TYPE##XMLConverter))); \
  

#define ADD_TEMPLATED_NUMBER_DEPS(T) \
  ADD_NUMBER_VISUAL_DEP(T); \
  ADD_RANGE_VALIDATOR_DEP(T); \
  ADD_NUMBER_ARRAY_LENGTH_DEP_GROUP(T);

#define ADD_NUMBER_VISUAL_DEP(T) \
  \
  masterMap.insert( \
    ConverterPair( \
      DummyObjectGetter<NumberVisualDependency< T > >:: \
      getDummyObject()->getTypeAttributeValue(), \
      rcp(new NumberVisualDependencyXMLConverter< T >)));

#define ADD_RANGE_VALIDATOR_DEP(T) \
  \
  masterMap.insert( \
    ConverterPair( \
      DummyObjectGetter<RangeValidatorDependency< T > >:: \
        getDummyObject()->getTypeAttributeValue(), \
      rcp(new RangeValidatorDependencyXMLConverter< T >)));

#define ADD_NUMBER_ARRAY_LENGTH_DEP(DEPENDEE_TYPE , DEPENDENT_TYPE) \
  masterMap.insert( \
    ConverterPair( \
      DummyObjectGetter<NumberArrayLengthDependency< \
        DEPENDEE_TYPE , DEPENDENT_TYPE > >::getDummyObject()->getTypeAttributeValue(), \
        rcp(new NumberArrayLengthDependencyXMLConverter< \
        DEPENDEE_TYPE , DEPENDENT_TYPE >)));

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
#define ADD_NUMBER_ARRAY_LENGTH_DEP_GROUP(DEPENDEE_TYPE) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , std::string) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , int) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , long long int) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , double) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , float) 
#else
#define ADD_NUMBER_ARRAY_LENGTH_DEP_GROUP(DEPENDEE_TYPE) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , std::string) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , int) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , double) \
  ADD_NUMBER_ARRAY_LENGTH_DEP( DEPENDEE_TYPE , float) 
#endif

  

void DependencyXMLConverterDB::addConverter(
  Dependency& dependency,
  RCP<DependencyXMLConverter> converterToAdd)
{
  getConverterMap().insert(
    ConverterPair(dependency.getTypeAttributeValue(), converterToAdd));
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
  if(masterMap.size() == 0){
    ADD_TEMPLATED_NUMBER_DEPS(int);
    ADD_NUMBER_VISUAL_DEP(float); 
    ADD_RANGE_VALIDATOR_DEP(float); 
    ADD_NUMBER_VISUAL_DEP(double); 
    ADD_RANGE_VALIDATOR_DEP(double); 

    #ifdef HAVE_TEUCHOS_LONG_LONG_INT
    ADD_TEMPLATED_NUMBER_DEPS(long long int);
    #endif // HAVE_TEUCHOS_LONG_LONG_INT

    ADD_DEP_TO_MAP(StringVisualDependency)
    ADD_DEP_TO_MAP(StringVisualDependency)
    ADD_DEP_TO_MAP(BoolValidatorDependency)
    ADD_DEP_TO_MAP(BoolVisualDependency)
    ADD_DEP_TO_MAP(ConditionVisualDependency)
  }
  return masterMap;
}


} // namespace Teuchos
