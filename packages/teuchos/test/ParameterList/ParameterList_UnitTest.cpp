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
	std::cout << "Are the written and read parameter lists the same? :";
	TEST_ASSERT(haveSameValues(myList, *readInPL));
}

TEUCHOS_UNIT_TEST(Teuchos_ParameterList, validatorXMLConverters){
	std::string xmlFileName = "ValidatorsConverterList.xml";
	std::string defaultParameterName = "default int";
	std::string minmaxParameterName = "min max int";
	Teuchos::ParameterList myList;
	Teuchos::RCP<EnhancedNumberValidator<int> > defaultIntValidator = rcp( new EnhancedNumberValidator<int>());
	Teuchos::RCP<EnhancedNumberValidator<int> > minMaxIntValidator = rcp( new EnhancedNumberValidator<int>(0,10));
	myList.set(defaultParameterName, 6, "int parameter with default validator", defaultIntValidator);
	myList.set(minmaxParameterName, 10, "int parameter with min and max validator", minMaxIntValidator);

	std::cout << "\n";
	std::cout << "XML Representation of Parameter List: \n";
	Teuchos::XMLParameterListWriter plWriter;
	XMLObject xmlPL = plWriter.toXML(myList);
	std::cout << xmlPL;

	std::cout << "writing xml to file...\n";
	writeParameterListToXmlFile(myList, xmlFileName);
	std::cout << "reading xml from file..\n";
	Teuchos::RCP<Teuchos::ParameterList> readInPL = getParametersFromXmlFile(xmlFileName);
	std::cout << "Are the written and read parameter lists the same? : ";
	TEST_ASSERT(haveSameValues(myList, *readInPL));
	std::cout << "Are the validators where they're suppoed to be? : \n";
	TEST_EQUALITY(
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator<int> >(readInPL->getEntry(defaultParameterName).validator())->getMin(),
		Teuchos::rcp_static_cast<const Teuchos::EnhancedNumberValidator<int> >(myList.getEntry(defaultParameterName).validator())->getMin()
	);



}

} // namespace Teuchos



