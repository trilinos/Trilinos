//@HEADER
// ***********************************************************************
//
//                           Teuchos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_as.hpp"


namespace Teuchos {

class UNDEFINED_PARAMETERENTRY_VALIDATOR : public ParameterEntryValidator
{
  
  public:

  void printDoc(const std::string& docString, std::ostream& out){}

  ValidStringsList validStringValues() const{
    return rcp(new tuple<std::string>(""));
  
  void validate(
    ParameterEntry  const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const {}

};


#define ADD_TYPE_PARAMETER(T,VALUE) \
  myList.set( #T , as<T>(VALUE));
const int g_arraySize = 5;
#define ADD_ARRAY_TYPE_PARAMETER(T,VALUE) \
  myList.set( #T " Array", Array< T >(g_arraySize, ( T ) VALUE ));
#define ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(T,VALUE) \
  ADD_TYPE_PARAMETER(T,VALUE); \
  ADD_ARRAY_TYPE_PARAMETER(T,VALUE);


/*
TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentNames ) {
  // Lists with different names should not be equal
  ParameterList A("Tom");
  ParameterList B("Bob");
  TEST_ASSERT( A != B );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( A != B );
}

TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentNames ) {
  ParameterList A("Julie");
  ParameterList B("Shannon");
  TEST_ASSERT( !haveSameValues(A,B) );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
}
*/


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityWithEmpty )
{
  // An empty list should not be equal to a full list
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( A == B );
  A.set("Hello","World");
  TEST_ASSERT( A != B );
  B.set("Hello","World");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentSublistNames )
{
  // Sublists with different names should not be equal
  ParameterList A;
  ParameterList B;
  A.sublist("Bob");
  B.sublist("Tom");
  TEST_ASSERT( A != B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentLengths )
{
  ParameterList A;
  ParameterList B;
  A.set("A","a");
  A.set("B","b");
  A.set("C","c");

  B.set("A","a");
  B.set("B","b");

  TEST_ASSERT( A != B );

  B.set("C","c");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesWithEmpty )
{
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( haveSameValues(A,B) );
  A.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
  B.set("Hello","World");
  TEST_ASSERT( haveSameValues(A,B) );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentSublistNames )
{
  ParameterList A;
  ParameterList B;
  A.sublist("Smith").set("People",4);
  B.sublist("Jones").set("People",4);
  TEST_ASSERT( !haveSameValues(A,B) ); // sublist names matter
}

TEUCHOS_UNIT_TEST(Teuchos_ParameterList, ADD_TYPE_AND_ARRAY_TYPE_PARAMETER)
{
  ParameterList myList;
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(int, 2);
  TEST_EQUALITY( getParameter<int>(myList, "int"), 2 );
  TEST_EQUALITY( getParameter<Array<int> >(myList, "int Array"),
    Array<int>(g_arraySize, as<int>(2)) );
}


TEUCHOS_UNIT_TEST(Teuchos_ParameterList, parameterEntryXMLConverters)
{
  std::string xmlFileName = "PEConvertersList.xml";
  ParameterList myList;
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(int, 2);
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(unsigned int, 3);
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(short int, 4);
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(unsigned short int, 5);
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(long int, 6);
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(unsigned long int, 7);
  #ifdef HAVE_TEUCHOS_LONG_LONG_INT
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(long long int, 8);
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(unsigned long long int, 9);
  #endif //HAVE_TEUCHOS_LONG_LONG_INT
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(double, 10.0);
  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(float, 11.0);

  ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(std::string, "hello");

  ADD_TYPE_PARAMETER(char, 'a');
  ADD_TYPE_PARAMETER(bool, true);

  out << "\nwriting xml to file...\n";
  writeParameterListToXmlFile(myList, xmlFileName);
  out << "reading xml from file...\n";
  RCP<ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
  TEST_ASSERT(haveSameValues(myList, *readInPL));
  if (!success) {
    out << "\nmyList:\n";
    myList.print(out);
    out << "\n*readInPL:\n";
    readInPL->print(out);
  }
}


TEUCHOS_UNIT_TEST(Teuchos_ParameterList, converterExceptions)
{

  TEST_THROW(RCP<ParameterList>
    missingValidatorList = getParametersFromXmlFile("MissingValidator.xml"),
    MissingValidatorDefinitionException);
 
  TEST_THROW(RCP<ParameterList>
    missingPrototypeList = getParametersFromXmlFile("MissingPrototypeValidator.xml"),
	MissingValidatorDefinitionException);

  TEST_THROW(RCP<ParameterList>
    badRootElementList = getParametersFromXmlFile("BadRootElement.xml"),
    BadXMLParameterListRootElementException);

  TEST_THROW(RCP<ParameterList>
    badParameterListElementList = getParametersFromXmlFile("BadParameterListElement.xml"),
    BadParameterListElementException);

  TEST_THROW(RCP<ParameterList>
    noNameAttributeList = getParametersFromXmlFile("NoNameAttribute.xml"),
    NoNameAttributeExecption);

  TEST_THROW(RCP<ParameterList>
    noTypeAttributeList = getParametersFromXmlFile("NoTypeAttribute.xml"),
    NoValueAttributeExecption);
  
  TEST_THROW(RCP<ParameterList>
    noValueAttributeList = getParametersFromXmlFile("NoValueAttribute.xml"),
    NoValueAttributeExecption);

  TEST_THROW(RCP<ParameterList>
    conflicitingValiIdsList = getParametersFromXmlFile("ConflictingValidatorIDs.xml"),
    DuplicateValidatorIDsException);

  TEST_THROW(RCP<ParameterList>
    stringValidatorBadTagList = getParametersFromXmlFile("StringValidatorBadTag.xml"),
    BadTagException);

  TEST_THROW(RCP<ParameterList>
    stringValidatorBadTagList = getParametersFromXmlFile("StringToIntegralValidatorBadTag.xml"),
    BadTagException);

  TEST_THROW(RCP<ParameterList>
    badParameterEntryConverterList = getParametersFromXmlFile("CantFindParameterEntryConverter.xml"),
	CantFindParameterEntryConverterException);

  UNDEFINED_PARAMETERENTRY_VALIDATOR badValidator;
  TEST_THROW(ValidatorXMLConverterDB::getConverter(badValidator), CantFindValidatorConverterException);

  #ifdef HAVE_TEUCHOS_DEBUG
  StandardTemplatedParameterConverter<int> intConverter;
  StandardTemplatedParameterConverter<float> floatConverter;
  ParameterEntry floatParameter("float", (float)3.0);
  TEST_THROW(intConverter.fromParameterEntryToXML(doubleParameter, "blah"), BadParameterEntryXMLConverterTypeExecption);

  XMLObject floatXML = floatConverter.fromParameterEntryToXML(floatParameter, "float");
  TEST_THROW(intConverter.fromXMLtoParameterEntry(floatXML), BadParameterEntryXMLConverterTypeExecption);

  IDtoValidatorMap dummyIDtoVMap;
  ValidatortoIDMap dummyVtoIDMap;
  StringValidatorXMLConverter stringConverter;
  AnyNumberValidatorXMLConverter anyNumberConverter;
  AnyNumberParameterEntryValidator anyNumberValidator;
  TEST_THROW(stringConverter.fromValidatortoXML(anyNumberValidator, dummyIDtoVMap), BadValidatorXMLConverterException);
  XMLObject anyNumberXML = anyNumberConverter(anyNumberValidator);
  TEST_THROW(stringConverter.fromXMLtoValidator(anyNumberXML, dummyVtoIDMap), BadValidatorXMLConverterException);
  #endif

}


TEUCHOS_UNIT_TEST(Teuchos_ParameterList, fileNameValidatorConverter)
{
  std::string xmlFileName = "FileNameValidatorList.xml";
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "non default";

  RCP<FileNameValidator> defaultValidator =
    rcp(new FileNameValidator);
  RCP<FileNameValidator> nonDefaultValidator =
    rcp(new FileNameValidator(true));
  ParameterList myList("FileName Validator List");
  myList.set("default", "blah.txt", "parameter for default validator",
    defaultValidator);
  myList.set("non default", "blah.txt", "parameter for non default validator",
    nonDefaultValidator);

  writeParameterListToXmlFile(myList, xmlFileName);
  RCP<ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);

  RCP<const FileNameValidator> readinDefault =
    rcp_static_cast<const FileNameValidator>(
      readInPL->getEntry(defaultParameterName).validator());
  TEST_EQUALITY(readinDefault->fileMustExist(), defaultValidator->fileMustExist());

  RCP<const FileNameValidator> readinNonDefault =
    rcp_static_cast<const FileNameValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator());
  TEST_EQUALITY(readinNonDefault->fileMustExist(), nonDefaultValidator->fileMustExist());
}


// 2010/07/30: rabartl: I am a little worried about all of the static casts.
// Can you not use dynamic casts?  You can only use (rcp_)dynamic_cast when
// the type has at least one virtual function (which all Validators should).
// 
// 2010/08/03 I'm using static casts because from what I understand
// they are faster than dynamic casts. there should be no danger because 
// I know these are valid casts.


TEUCHOS_UNIT_TEST(Teuchos_ParameterList, stringValidatorConverter)
{
  std::string xmlFileName = "StringValidatorList.xml";
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "non default";

  RCP<StringValidator> nonDefaultValidator = rcp(
    new StringValidator(tuple<std::string>("value1", "cheese", "kurtis", "is", "awesome")));
  ParameterList myList("String Validator List");
  myList.set("non default", "kurtis", "parameter for non default validator",
    nonDefaultValidator);

  writeParameterListToXmlFile(myList, xmlFileName);
  RCP<ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);

  RCP<const StringValidator> readinNonDefault =
    rcp_static_cast<const StringValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator());
  TEST_COMPARE_ARRAYS(*(readinNonDefault->validStringValues()),
    *(nonDefaultValidator->validStringValues()));
}


// 2010/07/30: rabartl: Could most of these unit tests be made to run faster
// by reading and writing XML from strings and not going to files?  Is there a
// need to go and actually write files?


TEUCHOS_UNIT_TEST(Teuchos_ParameterList, anynumberValidatorConverter)
{
  std::string xmlFileName = "AnyNumberValidatorList.xml";
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "preferred and accepted";
  RCP<AnyNumberParameterEntryValidator> defaultValidator =
    rcp(new AnyNumberParameterEntryValidator());
  AnyNumberParameterEntryValidator::AcceptedTypes acceptedTypes;
  acceptedTypes.allowDouble(false);
  RCP<AnyNumberParameterEntryValidator> nonDefaultValidator =
    rcp(
      new AnyNumberParameterEntryValidator(
        AnyNumberParameterEntryValidator::PREFER_INT,
        acceptedTypes
        )
      );

  ParameterList myList("AnyNumberValidatorList");
  myList.set(defaultParameterName, 10.0,
    "A parameter with the default AnyNumberValidator on it", defaultValidator);
  myList.set(nonDefaultParameterName, 1, 
    "A prameter with an AnyNumberValidator on it that has the preferred and accepted types differnet from the default",
    nonDefaultValidator);
  writeParameterListToXmlFile(myList, xmlFileName);
  RCP<ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
  
  RCP<const AnyNumberParameterEntryValidator> readinDefaultValidator =
    rcp_static_cast<const AnyNumberParameterEntryValidator>(
      readInPL->getEntry(defaultParameterName).validator());
  TEST_EQUALITY(readinDefaultValidator->isDoubleAllowed(),
    defaultValidator->isDoubleAllowed());
  TEST_EQUALITY(readinDefaultValidator->isIntAllowed(),
    defaultValidator->isIntAllowed());
  TEST_EQUALITY(readinDefaultValidator->isStringAllowed(),
    defaultValidator->isStringAllowed());
  TEST_EQUALITY(readinDefaultValidator->getPreferredType(),
    defaultValidator->getPreferredType());

  RCP<const AnyNumberParameterEntryValidator> readinNonDefaultValidator =
    rcp_static_cast<const AnyNumberParameterEntryValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator());
  TEST_EQUALITY(readinNonDefaultValidator->isDoubleAllowed(),
    nonDefaultValidator->isDoubleAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->isIntAllowed(),
    nonDefaultValidator->isIntAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->isStringAllowed(),
    nonDefaultValidator->isStringAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->getPreferredType(),
    nonDefaultValidator->getPreferredType());
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_ParameterList, EnhancedNumberValidatorConverter, T)
{
  std::string xmlFileName = TypeNameTraits<T>::name() + "EnhancedValidatorList.xml";
  std::string defaultParameterName = "default";
  std::string minmaxParameterName = "min max";
  ParameterList myList;
  RCP<EnhancedNumberValidator< T > > defaultValidator =
    rcp( new EnhancedNumberValidator< T >());
  RCP<EnhancedNumberValidator< T > > minMaxValidator =
    rcp( new EnhancedNumberValidator< T >(0,10));
  myList.set(defaultParameterName, ( T )6, "parameter with default validator",
    defaultValidator);
  myList.set(minmaxParameterName, ( T )10, "parameter with min and max validator",
    minMaxValidator);
  writeParameterListToXmlFile(myList, xmlFileName);
  RCP<ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator())->getMin(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator())->getMin()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator())->getMax(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator())->getMax()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator())->getStep(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator())->getStep()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator())->getPrecision(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator())->getPrecision()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator())->hasMin(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator())->hasMin()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator())->hasMax(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator())->hasMax()
  );

  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator())->getMin(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator())->getMin()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator())->getMax(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator())->getMax()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator())->getStep(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator())->getStep()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator())->getPrecision(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator())->getPrecision()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator())->hasMin(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator())->hasMin()
  );
  TEST_EQUALITY(
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator())->hasMax(),
    rcp_static_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator())->hasMax()
  );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_ParameterList, NumberArrayValidatorConverterTest, T)
{
  std::string xmlFileName = TypeNameTraits<T>::name() + "ArrayConverterList.xml";
  std::string arrayParameterName = "array";
  ParameterList myList;

  const T arrayValidatorLen = as<T>(11);
  RCP<ArrayNumberValidator< T > > arrayValidator =
    rcp(new ArrayNumberValidator< T >(
          rcp(new EnhancedNumberValidator<T>(as<T>(0), arrayValidatorLen))));
  myList.set(arrayParameterName,   Array< T >(4, 10), "array parameter", arrayValidator);
  writeParameterListToXmlFile(myList, xmlFileName);

  RCP<ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
  RCP<const EnhancedNumberValidator< T > > readInPrototypeValidator =
    rcp_static_cast<const ArrayNumberValidator< T > >(
      readInPL->getEntry(arrayParameterName).validator())->getPrototype();
  RCP<const EnhancedNumberValidator< T > > actualPrototypeValidator =
    rcp_static_cast<const ArrayNumberValidator< T > >(
      myList.getEntry(arrayParameterName).validator())->getPrototype();

  TEST_EQUALITY(
    readInPrototypeValidator->getMin(),
    actualPrototypeValidator->getMin()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->getMax(),
    actualPrototypeValidator->getMax()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->getStep(),
    actualPrototypeValidator->getStep()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->getPrecision(),
    actualPrototypeValidator->getPrecision()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->hasMin(),
    actualPrototypeValidator->hasMin()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->hasMax(),
    actualPrototypeValidator->hasMax()
  );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_ParameterList, StringToIntegralConverterTest, T)
{
  std::string xmlFileName = TypeNameTraits<T>::name() + "STIConverterList.xml";
  std::string defaultStringToIntegralParameterName = "defaultsti";
  std::string stringToIntegralParameterName = "sti";
  ParameterList myList;
  RCP<StringToIntegralParameterEntryValidator< T > > defaultStiValidator = rcp(
    new StringToIntegralParameterEntryValidator< T >(
      tuple<std::string>("value1", "value2", "value3"), stringToIntegralParameterName));
  RCP<StringToIntegralParameterEntryValidator< T > > stiValidator = rcp(
    new StringToIntegralParameterEntryValidator< T >(
      tuple<std::string>("value3", "value4", "value5"), 
      tuple<std::string>("the third value", "the fourth value", "the fifth value"),
      tuple< T >(3,4,5),
      stringToIntegralParameterName));
  myList.set(defaultStringToIntegralParameterName,
    "value1", "parameter with default sti validator", defaultStiValidator);
  myList.set(stringToIntegralParameterName, "value3", "parameter with sti validator",
    stiValidator);

  writeParameterListToXmlFile(myList, xmlFileName);
  RCP<ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);


  RCP<const StringToIntegralParameterEntryValidator< T > > readInDefaultStiValidator =
    rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(
      readInPL->getEntry(defaultStringToIntegralParameterName).validator());
  RCP<const StringToIntegralParameterEntryValidator< T > > readInStiValidator =
    rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(
      readInPL->getEntry(stringToIntegralParameterName).validator());
  RCP<const StringToIntegralParameterEntryValidator< T > > actualDefaultStiValidator =
    rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(
      myList.getEntry(defaultStringToIntegralParameterName).validator());
  RCP<const StringToIntegralParameterEntryValidator< T > > actualStiValidator =
    rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(
      myList.getEntry(stringToIntegralParameterName).validator());

  Array<std::string> readInDefaultValidStrings =
    *(readInDefaultStiValidator->validStringValues());
  Array<std::string> actualDefaultValidStrings =
    *(actualDefaultStiValidator->validStringValues());
  TEST_COMPARE_ARRAYS(readInDefaultValidStrings, actualDefaultValidStrings);

  TEST_ASSERT(readInDefaultStiValidator->getStringDocs().is_null());
  TEST_EQUALITY( readInDefaultStiValidator->getDefaultParameterName(),
    actualDefaultStiValidator->getDefaultParameterName());
  for(int i=0; i<actualDefaultValidStrings.size(); ++i){
    TEST_EQUALITY(actualDefaultStiValidator->getIntegralValue(actualDefaultValidStrings[i]),
      readInDefaultStiValidator->getIntegralValue(actualDefaultValidStrings[i]));
  }

  Array<std::string> readInValidStrings = *(readInStiValidator->validStringValues());
  Array<std::string> actualValidStrings = *(actualStiValidator->validStringValues());
  TEST_COMPARE_ARRAYS(readInValidStrings, actualValidStrings);

  TEST_COMPARE_ARRAYS(*(readInStiValidator->getStringDocs()),
    *(actualStiValidator->getStringDocs()));
  TEST_EQUALITY( readInStiValidator->getDefaultParameterName(),
    actualStiValidator->getDefaultParameterName());
  for(int i=0; i<actualValidStrings.size(); ++i){
    TEST_EQUALITY(actualStiValidator->getIntegralValue(actualValidStrings[i]),
      readInStiValidator->getIntegralValue(actualValidStrings[i]));
  }

}

TEUCHOS_UNIT_TEST(Teuchos_ParameterList, existingPrototypeTest){
  ParameterList pl("ExsitingPrototypeList");
  RCP<StringValidator> stringVali = rcp(new StringValidator());
  RCP<ArrayValidator<StringValidator, std::string> > arrayStringVali 
    = rcp(new ArrayValidator<StringValidator, std::string>(stringVali));
  Array<std::string> strArray = tuple<std::string>("blah", "blah", "blah");
  pl.set("string param", "hi", "a string param", stringVali);
  pl.set("string array param", strArray, "a string array parameter", arrayStringVali);
  XMLParameterListWriter plWriter;
  XMLObject plXMLOut = plWriter.toXML(pl);
  std::string plStrOut = toString(plXMLOut);
  RCP<StringInputStream> plInStream = rcp( new StringInputStream(plStrOut));
  XMLParser xmlParser(plInStream);
  XMLObject plXMLIn = xmlParser.parse();
  XMLParameterListReader plReader;
  ParameterList plIn = plReader.toParameterList(plXMLIn);
  RCP<ArrayValidator<StringValidator, std::string> > inArrayValidator = rpc_dynamic_cast(plIn.getEntry("string array param").validator());
  TEST_ASSERT(*(plIn.getEntry().validator()) == *(inArrayValidator->getPrototype()));
}


#define FULL_NUMBER_TYPE_TEST( T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_ParameterList, EnhancedNumberValidatorConverter, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_ParameterList, NumberArrayValidatorConverterTest, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_ParameterList, StringToIntegralConverterTest, T ) 


typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;


FULL_NUMBER_TYPE_TEST(int)
FULL_NUMBER_TYPE_TEST(uint)
FULL_NUMBER_TYPE_TEST(short)
FULL_NUMBER_TYPE_TEST(ushort)
FULL_NUMBER_TYPE_TEST(long)
FULL_NUMBER_TYPE_TEST(ulong)
FULL_NUMBER_TYPE_TEST(float)
FULL_NUMBER_TYPE_TEST(double)
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int llint;
typedef unsigned long long int ullint;
FULL_NUMBER_TYPE_TEST(llint)
FULL_NUMBER_TYPE_TEST(ullint)
#endif


} // namespace Teuchos



