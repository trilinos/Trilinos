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

namespace Teuchos {

#define ADD_TYPE_PARAMETER(T,VALUE) myList.set( #T , ( T ) VALUE);
#define ADD_ARRAY_TYPE_PARAMETER(T,VALUE)  myList.set( #T " Array", Array< T >(5, ( T ) VALUE ));
#define ADD_TYPE_AND_ARRAY_TYPE_PARAMETER(T,VALUE) \
	ADD_TYPE_PARAMETER(T,VALUE); \
	ADD_ARRAY_TYPE_PARAMETER(T,VALUE);

/*
TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentNames ) {
  // Lists with different names should not be equal
  Teuchos::ParameterList A("Tom");
  Teuchos::ParameterList B("Bob");
  TEST_ASSERT( A != B );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( A != B );
}

TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentNames ) {
  Teuchos::ParameterList A("Julie");
  Teuchos::ParameterList B("Shannon");
  TEST_ASSERT( !haveSameValues(A,B) );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
}
*/

TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityWithEmpty ) {
  // An empty list should not be equal to a full list
  Teuchos::ParameterList A;
  Teuchos::ParameterList B;
  TEST_ASSERT( A == B );
  A.set("Hello","World");
  TEST_ASSERT( A != B );
  B.set("Hello","World");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentSublistNames ) {
  // Sublists with different names should not be equal
  Teuchos::ParameterList A;
  Teuchos::ParameterList B;
  A.sublist("Bob");
  B.sublist("Tom");
  TEST_ASSERT( A != B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentLengths ) {
  Teuchos::ParameterList A;
  Teuchos::ParameterList B;
  A.set("A","a");
  A.set("B","b");
  A.set("C","c");

  B.set("A","a");
  B.set("B","b");

  TEST_ASSERT( A != B );

  B.set("C","c");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesWithEmpty ) {
  Teuchos::ParameterList A;
  Teuchos::ParameterList B;
  TEST_ASSERT( haveSameValues(A,B) );
  A.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
  B.set("Hello","World");
  TEST_ASSERT( haveSameValues(A,B) );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentSublistNames ) {
  Teuchos::ParameterList A;
  Teuchos::ParameterList B;
  A.sublist("Smith").set("People",4);
  B.sublist("Jones").set("People",4);
  TEST_ASSERT( !haveSameValues(A,B) ); // sublist names matter
}

TEUCHOS_UNIT_TEST(Teuchos_ParameterList, parameterEntryXMLConverters){
	std::string xmlFileName = "PEConvertersList.xml";
	Teuchos::ParameterList myList;
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

	std::cout << "\n";
	std::cout << "XML Representation of Parameter List: \n";
	Teuchos::XMLParameterListWriter plWriter;
	XMLObject xmlPL = plWriter.toXML(myList);
	std::cout << xmlPL;

	std::cout << "writing xml to file...\n";
	writeParameterListToXmlFile(myList, xmlFileName);
	std::cout << "reading xml from file...\n";
	Teuchos::RCP<Teuchos::ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
	TEST_ASSERT(haveSameValues(myList, *readInPL));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_ParameterList, EnhancedNumberValidatorConverter, T){
	std::string xmlFileName = TypeNameTraits<T>::name() + "EnhancedValidatorList.xml";
	std::string defaultParameterName = "default";
	std::string minmaxParameterName = "min max";
	Teuchos::ParameterList myList;
	Teuchos::RCP<EnhancedNumberValidator< T > > defaultValidator = rcp( new EnhancedNumberValidator< T >());
	Teuchos::RCP<EnhancedNumberValidator< T > > minMaxValidator = rcp( new EnhancedNumberValidator< T >(0,10));
	myList.set(defaultParameterName, ( T )6, "parameter with default validator", defaultValidator);
	myList.set(minmaxParameterName, ( T )10, "parameter with min and max validator", minMaxValidator);
	writeParameterListToXmlFile(myList, xmlFileName);
	Teuchos::RCP<Teuchos::ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(defaultParameterName).validator())->getMin(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(defaultParameterName).validator())->getMin()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(defaultParameterName).validator())->getMax(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(defaultParameterName).validator())->getMax()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(defaultParameterName).validator())->getStep(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(defaultParameterName).validator())->getStep()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(defaultParameterName).validator())->getPrecision(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(defaultParameterName).validator())->getPrecision()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(defaultParameterName).validator())->hasMin(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(defaultParameterName).validator())->hasMin()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(defaultParameterName).validator())->hasMax(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(defaultParameterName).validator())->hasMax()
	);

	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(minmaxParameterName).validator())->getMin(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(minmaxParameterName).validator())->getMin()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(minmaxParameterName).validator())->getMax(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(minmaxParameterName).validator())->getMax()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(minmaxParameterName).validator())->getStep(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(minmaxParameterName).validator())->getStep()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(minmaxParameterName).validator())->getPrecision(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(minmaxParameterName).validator())->getPrecision()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(minmaxParameterName).validator())->hasMin(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(minmaxParameterName).validator())->hasMin()
	);
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(readInPL->getEntry(minmaxParameterName).validator())->hasMax(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator< T > >(myList.getEntry(minmaxParameterName).validator())->hasMax()
	);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_ParameterList, NumberArrayValidatorConverterTest, T){
	std::string xmlFileName = TypeNameTraits<T>::name() + "ArrayConverterList.xml";
	std::string arrayParameterName = "array";
	Teuchos::ParameterList myList;
	Teuchos::RCP<ArrayNumberValidator< T > > arrayValidator = rcp(new ArrayNumberValidator< T >(rcp(new EnhancedNumberValidator< T >(( T )0,( T )11))));
	myList.set(arrayParameterName, 	Teuchos::Array< T >(4, 10), "array parameter", arrayValidator);
	writeParameterListToXmlFile(myList, xmlFileName);
	Teuchos::RCP<Teuchos::ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
	Teuchos::RCP<const EnhancedNumberValidator< T > > readInPrototypeValidator = Teuchos::rcp_static_cast<const ArrayNumberValidator< T > >(readInPL->getEntry(arrayParameterName).validator())->getPrototype();
	Teuchos::RCP<const EnhancedNumberValidator< T > > actualPrototypeValidator = Teuchos::rcp_static_cast<const ArrayNumberValidator< T > >(myList.getEntry(arrayParameterName).validator())->getPrototype();
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

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_ParameterList, StringToIntegralConverterTest, T){
	std::string xmlFileName = TypeNameTraits<T>::name() + "STIConverterList.xml";
	std::string defaultStringToIntegralParameterName = "defaultsti";
	std::string stringToIntegralParameterName = "sti";
	Teuchos::ParameterList myList;
	Teuchos::RCP<StringToIntegralParameterEntryValidator< T > > defaultStiValidator = rcp(
		new StringToIntegralParameterEntryValidator< T >(tuple<std::string>("value1", "value2", "value3"), stringToIntegralParameterName));
	Teuchos::RCP<StringToIntegralParameterEntryValidator< T > > stiValidator = rcp(
		new StringToIntegralParameterEntryValidator< T >(
			tuple<std::string>("value3", "value4", "value5"), 
			tuple<std::string>("the third value", "the fourth value", "the fifth value"),
			tuple< T >(3,4,5),
			stringToIntegralParameterName));
	myList.set(defaultStringToIntegralParameterName,"value1", "parameter with default sti validator", defaultStiValidator);
	myList.set(stringToIntegralParameterName,"value3", "parameter with sti validator", stiValidator);

	writeParameterListToXmlFile(myList, xmlFileName);
	Teuchos::RCP<Teuchos::ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);


	Teuchos::RCP<const StringToIntegralParameterEntryValidator< T > > readInDefaultStiValidator = Teuchos::rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(readInPL->getEntry(defaultStringToIntegralParameterName).validator());
	Teuchos::RCP<const StringToIntegralParameterEntryValidator< T > > readInStiValidator = Teuchos::rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(readInPL->getEntry(stringToIntegralParameterName).validator());
	Teuchos::RCP<const StringToIntegralParameterEntryValidator< T > > actualDefaultStiValidator = Teuchos::rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(myList.getEntry(defaultStringToIntegralParameterName).validator());
	Teuchos::RCP<const StringToIntegralParameterEntryValidator< T > > actualStiValidator = Teuchos::rcp_static_cast<const StringToIntegralParameterEntryValidator< T > >(myList.getEntry(stringToIntegralParameterName).validator());

	Teuchos::Array<std::string> readInDefaultValidStrings = *(readInDefaultStiValidator->validStringValues());
	Teuchos::Array<std::string> actualDefaultValidStrings = *(actualDefaultStiValidator->validStringValues());
	TEST_COMPARE_ARRAYS(readInDefaultValidStrings, actualDefaultValidStrings);

	TEST_ASSERT(readInDefaultStiValidator->getStringDocs().is_null());
	TEST_EQUALITY( readInDefaultStiValidator->getDefaultParameterName(), actualDefaultStiValidator->getDefaultParameterName());
	for(int i=0; i<actualDefaultValidStrings.size(); ++i){
		TEST_EQUALITY(actualDefaultStiValidator->getIntegralValue(actualDefaultValidStrings[i]), readInDefaultStiValidator->getIntegralValue(actualDefaultValidStrings[i]));
	}

	Teuchos::Array<std::string> readInValidStrings = *(readInStiValidator->validStringValues());
	Teuchos::Array<std::string> actualValidStrings = *(actualStiValidator->validStringValues());
	TEST_COMPARE_ARRAYS(readInValidStrings, actualValidStrings);

	TEST_COMPARE_ARRAYS(*(readInStiValidator->getStringDocs()), *(actualStiValidator->getStringDocs()));
	TEST_EQUALITY( readInStiValidator->getDefaultParameterName(), actualStiValidator->getDefaultParameterName());
	for(int i=0; i<actualValidStrings.size(); ++i){
		TEST_EQUALITY(actualStiValidator->getIntegralValue(actualValidStrings[i]), readInStiValidator->getIntegralValue(actualValidStrings[i]));
	}

}

#define FULL_NUMBER_TYEP_TEST( T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_ParameterList, EnhancedNumberValidatorConverter, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_ParameterList, NumberArrayValidatorConverterTest, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_ParameterList, StringToIntegralConverterTest, T ) 

FULL_NUMBER_TYEP_TEST(int)
FULL_NUMBER_TYEP_TEST(short)
FULL_NUMBER_TYEP_TEST(long)
FULL_NUMBER_TYEP_TEST(float)
FULL_NUMBER_TYEP_TEST(double)

} // namespace Teuchos



