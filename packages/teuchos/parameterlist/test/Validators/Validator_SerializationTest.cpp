// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardValidatorXMLConverters.hpp"

#include "Teuchos_XMLParameterListTestHelpers.hpp"

using Teuchos::writeThenReadPL;


namespace Teuchos {

class UNDEFINED_PARAMETERENTRY_VALIDATOR : public ParameterEntryValidator
{

  public:

  void printDoc(const std::string& docString, std::ostream& out) const {}

  ValidStringsList validStringValues() const{
    return rcp(new Array<std::string>(1,""));
  }

  void validate(
    ParameterEntry  const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const {}

  const std::string getXMLTypeName() const{
    return "UNDEFINEDTYPE";
  }

};

TEUCHOS_UNIT_TEST(Teuchos_Validator, exceptionTests)
{
  ValidatorXMLConverterDB::printKnownConverters(out);
  out << std::endl;

  UNDEFINED_PARAMETERENTRY_VALIDATOR badValidator;
  TEST_THROW(ValidatorXMLConverterDB::getConverter(badValidator), CantFindValidatorConverterException);

  TEST_THROW(RCP<ParameterList>
    missingValidatorList = getParametersFromXmlFile("MissingValidator.xml"),
    MissingValidatorDefinitionException);

  TEST_THROW(RCP<ParameterList>
    missingPrototypeList = getParametersFromXmlFile("MissingPrototypeValidator.xml"),
        MissingValidatorDefinitionException);

  TEST_THROW(RCP<ParameterList>
    conflicitingValiIdsList = getParametersFromXmlFile("ConflictingValidatorIDs.xml"),
    DuplicateValidatorIDsException);

  TEST_THROW(RCP<ParameterList>
    stringValidatorBadTagList = getParametersFromXmlFile("StringValidatorBadTag.xml"),
    BadTagException);

  TEST_THROW(RCP<ParameterList>
    stringValidatorBadTagList = getParametersFromXmlFile("StringToIntegralValidatorBadTag.xml"),
    BadTagException);

  #ifdef HAVE_TEUCHOS_DEBUG

  StringValidatorXMLConverter stringConverter;
  AnyNumberValidatorXMLConverter anyNumberConverter;
  ValidatortoIDMap writerDummyMap;
  IDtoValidatorMap readerDummyMap;
  RCP<AnyNumberParameterEntryValidator> anyNumberValidator =
    anyNumberParameterEntryValidator();
  writerDummyMap.insert(anyNumberValidator);
  TEST_THROW(
    stringConverter.fromValidatortoXML(anyNumberValidator, writerDummyMap),
    BadValidatorXMLConverterException);
  XMLObject anyNumberXML =
    anyNumberConverter.fromValidatortoXML(anyNumberValidator, writerDummyMap);
  TEST_THROW(
    stringConverter.fromXMLtoValidator(anyNumberXML, readerDummyMap),
    BadValidatorXMLConverterException);

  #endif

}

TEUCHOS_UNIT_TEST(Teuchos_Validator, fileNameValidatorConverter)
{
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "non default";

  RCP<FileNameValidator> defaultValidator =
    rcp(new FileNameValidator);
  RCP<FileNameValidator> nonDefaultValidator =
    rcp(new FileNameValidator(true));
  ParameterList myList("FileName Validator List");
  myList.set("default", "", "parameter for default validator",
    defaultValidator);
  myList.set("non default", "blah.txt", "parameter for non default validator",
    nonDefaultValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const FileNameValidator> readinDefault =
    rcp_dynamic_cast<const FileNameValidator>(
      readInPL->getEntry(defaultParameterName).validator(), true);
  TEST_EQUALITY(readinDefault->fileMustExist(), defaultValidator->fileMustExist());

  RCP<const FileNameValidator> readinNonDefault =
    rcp_dynamic_cast<const FileNameValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator(), true);
  TEST_EQUALITY(readinNonDefault->fileMustExist(), nonDefaultValidator->fileMustExist());
}


TEUCHOS_UNIT_TEST(Teuchos_Validator, stringValidatorConverter)
{
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "non default";

  RCP<StringValidator> nonDefaultValidator = rcp(
    new StringValidator(tuple<std::string>("value1", "cheese", "kurtis", "is", "awesome")));
  ParameterList myList("String Validator List");
  myList.set("non default", "kurtis", "parameter for non default validator",
    nonDefaultValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const StringValidator> readinNonDefault =
    rcp_dynamic_cast<const StringValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator(), true);
  TEST_COMPARE_ARRAYS(*(readinNonDefault->validStringValues()),
    *(nonDefaultValidator->validStringValues()));
}


TEUCHOS_UNIT_TEST(Teuchos_Validator, boolValidatorConverter)
{
  std::string xmlFileName = "BoolValidatorList.xml";
  std::string boolParameterName = "boolParameterName";
  RCP<BoolParameterEntryValidator> boolValidator =
    rcp(new BoolParameterEntryValidator());

  ParameterList myList("BoolValidatorList");
  myList.set(boolParameterName, false,
    "A parameter with a BoolParameterEntryValidator validator.",
    boolValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const BoolParameterEntryValidator> readInBoolValidator =
    rcp_dynamic_cast<const BoolParameterEntryValidator>(
      readInPL->getEntry(boolParameterName).validator(), true);

  // to do - check any stuff we want to check
  // right now it doesn't have any settings
}

TEUCHOS_UNIT_TEST(Teuchos_Validator, anynumberValidatorConverter)
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

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const AnyNumberParameterEntryValidator> readinDefaultValidator =
    rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(
      readInPL->getEntry(defaultParameterName).validator(), true);
  TEST_EQUALITY(readinDefaultValidator->isDoubleAllowed(),
    defaultValidator->isDoubleAllowed());
  TEST_EQUALITY(readinDefaultValidator->isIntAllowed(),
    defaultValidator->isIntAllowed());
  TEST_EQUALITY(readinDefaultValidator->isStringAllowed(),
    defaultValidator->isStringAllowed());
  TEST_EQUALITY(readinDefaultValidator->getPreferredType(),
    defaultValidator->getPreferredType());

  RCP<const AnyNumberParameterEntryValidator> readinNonDefaultValidator =
    rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator(), true);
  TEST_EQUALITY(readinNonDefaultValidator->isDoubleAllowed(),
    nonDefaultValidator->isDoubleAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->isIntAllowed(),
    nonDefaultValidator->isIntAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->isStringAllowed(),
    nonDefaultValidator->isStringAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->getPreferredType(),
    nonDefaultValidator->getPreferredType());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Validator, EnhancedNumberValidatorConverter, T)
{
  std::string xmlFileName = TypeNameTraits<T>::name() + "EnhancedValidatorList.xml";
  std::string defaultParameterName = "default";
  std::string minmaxParameterName = "min max";
  std::string stepPrecParameterName = "step and prec";
  ParameterList myList;
  RCP<EnhancedNumberValidator< T > > defaultValidator =
    rcp( new EnhancedNumberValidator< T >());
  RCP<EnhancedNumberValidator< T > > minMaxValidator =
    rcp( new EnhancedNumberValidator< T >(0,10));
  RCP<EnhancedNumberValidator< T > > stepAndPrecValidator =
    rcp( new EnhancedNumberValidator< T >(0,10,4,4));
  myList.set(defaultParameterName, ( T )6, "parameter with default validator",
    defaultValidator);
  myList.set(minmaxParameterName, ( T )10, "parameter with min and max validator",
    minMaxValidator);
  myList.set(stepPrecParameterName, ( T )10, "parameter with min, max, "
    "step, and prec validator",
    stepAndPrecValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->getMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->getMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->getMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->getMax()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->getStep()
    ,
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->getStep()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(
        defaultParameterName).validator(), true)->getPrecision(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(
        defaultParameterName).validator(), true)->getPrecision()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->hasMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->hasMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->hasMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->hasMax()
  );

  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->getMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->getMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->getMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->getMax()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->getStep(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->getStep()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(
        minmaxParameterName).validator(), true)->getPrecision(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(
        minmaxParameterName).validator(), true)->getPrecision()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->hasMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->hasMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->hasMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->hasMax()
  );

  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(stepPrecParameterName).validator(), true)->getMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(stepPrecParameterName).validator(), true)->getMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(stepPrecParameterName).validator(), true)->getMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(stepPrecParameterName).validator(), true)->getMax()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(stepPrecParameterName).validator(), true)->getStep()
    ,
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(stepPrecParameterName).validator(), true)->getStep()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(
        stepPrecParameterName).validator(), true)->getPrecision(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(
        stepPrecParameterName).validator(), true)->getPrecision());

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Validator, NumberArrayValidatorConverterTest, T)
{
  std::string arrayParameterName = "array";
  ParameterList myList;

  const T arrayValidatorLen = as<T>(11);
  RCP<ArrayNumberValidator< T > > arrayValidator =
    rcp(new ArrayNumberValidator< T >(
      rcp(new EnhancedNumberValidator<T>(as<T>(0), arrayValidatorLen))));
  myList.set(arrayParameterName,
    Array< T >(4, 10), "array parameter", arrayValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const EnhancedNumberValidator< T > > readInPrototypeValidator =
    rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<T>, T > >(
      readInPL->getEntry(
        arrayParameterName).validator(), true)->getPrototype();
  RCP<const EnhancedNumberValidator< T > > actualPrototypeValidator =
    arrayValidator->getPrototype();

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

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Validator, TwoDArrayNumberValidatorConverterTest, T)
{
  std::string arrayParameterName = "array";
  ParameterList myList;

  const T arrayValidatorLen = as<T>(11);
  RCP<TwoDArrayNumberValidator< T > > arrayValidator =
    rcp(new TwoDArrayNumberValidator< T >(
      rcp(new EnhancedNumberValidator<T>(as<T>(0), arrayValidatorLen))));
  myList.set(arrayParameterName,
    TwoDArray< T >(4,4, 10), "array parameter", arrayValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const EnhancedNumberValidator< T > > readInPrototypeValidator =
    rcp_dynamic_cast<const TwoDArrayValidator<EnhancedNumberValidator<T>, T > >(
      readInPL->getEntry(
        arrayParameterName).validator(), true)->getPrototype();
  RCP<const EnhancedNumberValidator< T > > actualPrototypeValidator =
    arrayValidator->getPrototype();

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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Validator, StringToIntegralConverterTest, T)
{
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

  RCP<ParameterList> readInPL = writeThenReadPL(myList);


  RCP<const StringToIntegralParameterEntryValidator< T > >
  readInDefaultStiValidator =
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator< T > >(
      readInPL->getEntry(
        defaultStringToIntegralParameterName).validator(), true);
  RCP<const StringToIntegralParameterEntryValidator< T > >
  readInStiValidator =
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator< T > >(
      readInPL->getEntry(
        stringToIntegralParameterName).validator(), true);

  Array<std::string> readInDefaultValidStrings =
    *(readInDefaultStiValidator->validStringValues());
  Array<std::string> defaultValidStrings =
    *(defaultStiValidator->validStringValues());
  TEST_COMPARE_ARRAYS(readInDefaultValidStrings, defaultValidStrings);

  TEST_ASSERT(readInDefaultStiValidator->getStringDocs().is_null());
  TEST_EQUALITY( readInDefaultStiValidator->getDefaultParameterName(),
    defaultStiValidator->getDefaultParameterName());
  for(int i=0; i<defaultValidStrings.size(); ++i){
    TEST_EQUALITY(defaultStiValidator->getIntegralValue(defaultValidStrings[i]),
      readInDefaultStiValidator->getIntegralValue(defaultValidStrings[i]));
  }

  Array<std::string> readInValidStrings = *(readInStiValidator->validStringValues());
  Array<std::string> validStrings = *(stiValidator->validStringValues());
  TEST_COMPARE_ARRAYS(readInValidStrings, validStrings);

  TEST_COMPARE_ARRAYS(*(readInStiValidator->getStringDocs()),
    *(stiValidator->getStringDocs()));
  TEST_EQUALITY( readInStiValidator->getDefaultParameterName(),
    stiValidator->getDefaultParameterName());
  for(int i=0; i<validStrings.size(); ++i){
    TEST_EQUALITY(stiValidator->getIntegralValue(validStrings[i]),
      readInStiValidator->getIntegralValue(validStrings[i]));
  }

}

TEUCHOS_UNIT_TEST(Teuchos_Validator, existingPrototypeTest){
  ParameterList pl("ExsitingPrototypeList");
  RCP<StringValidator> stringVali = rcp(new StringValidator());
  RCP<ArrayValidator<StringValidator, std::string> > arrayStringVali
    = rcp(new ArrayValidator<StringValidator, std::string>(stringVali));
  Array<std::string> strArray = tuple<std::string>("blah", "blah", "blah");
  pl.set("string param", "hi", "a string param", stringVali);
  pl.set("string array param", strArray,
    "a string array parameter", arrayStringVali);
  RCP<ParameterList> readInPL = writeThenReadPL(pl);
  RCP<const ArrayValidator<StringValidator, std::string> >
    inArrayValidator =
    rcp_dynamic_cast<const ArrayValidator<StringValidator, std::string> >(
      readInPL->getEntry("string array param").validator(), true);
  TEST_ASSERT(readInPL->getEntry("string param").validator()
    == inArrayValidator->getPrototype());
}


#define FULL_NUMBER_TYPE_TEST( T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, EnhancedNumberValidatorConverter, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, NumberArrayValidatorConverterTest, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, StringToIntegralConverterTest, T )

#define NONINTEGRAL_NUMBER_TYPE_TEST( T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, EnhancedNumberValidatorConverter, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, NumberArrayValidatorConverterTest, T )

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;


FULL_NUMBER_TYPE_TEST(int)
NONINTEGRAL_NUMBER_TYPE_TEST(double)
NONINTEGRAL_NUMBER_TYPE_TEST(float)
typedef long long int llint;
FULL_NUMBER_TYPE_TEST(llint)

} // namespace Teuchos

