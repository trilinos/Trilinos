/*
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
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_as.hpp"
#include <iostream>

//ignore this for now
class CustomDataType{
public:
  CustomDataType():_theInt(0), _theString(""){}
  CustomDataType(int theInt, std::string theString):_theInt(theInt), _theString(theString){}

  void setInt(int theInt){
    _theInt = theInt;
  }

  void setString(std::string theString){
    _theString = theString;
  }

  int getInt() const{
    return _theInt;
  }

  std::string getString() const{
    return _theString;
  }

  bool operator==(const CustomDataType &other) const{
    return _theInt == other._theInt && _theString == other._theString;
  }

private:
  int _theInt;
  std::string _theString;
};

std::ostream& operator<<(std::ostream &out, const CustomDataType &object){
  out << object.getInt() << " " << object.getString();
  return out;
}

std::istream& operator>>(std::istream &in, CustomDataType &object){
  int theInt;
  std::string theString;
  in >> theInt;
  in >> theString;
  object.setInt(theInt);
  object.setString(theString);
  return in;
}

/**
 * this examples gives a demonstration of all the supporeted data types
 * for a parameter list. By supported we mean that they be serailized to
 * and from XML. You can use any data type you want as a parameter, but only
 * these types will be serailized to and from XML properly. In order
 * to serialize another data type, you must create and converter, something
 * which we will demonstrate in this example.
 */
int main(int argc, char* argv[])
{

  using Teuchos::tuple;
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  ParameterList myPL;

  //Basic data types
  myPL.set<int>("my int", 1);
  myPL.set<unsigned int>("my unsigned int", 1);
  myPL.set<short int>("my short int", 1);
  myPL.set<short>("my short", 1);
  myPL.set<unsigned short int>("my unsigned short int", 1);
  myPL.set<unsigned short>("my unsigned short", 1);
  myPL.set<long int>("my long int", 1);
  myPL.set<long>("my long", 1);
  myPL.set<unsigned long int>("my unsigned long int", 1);
  myPL.set<unsigned long>("my unsigned long", 1);
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  myPL.set<long long int>("my long long int", 1);
  myPL.set<long long>("my long long", 1);
  myPL.set<unsigned long long int>("my unsigned long long int", 1);
  myPL.set<unsigned long long>("my unsigned long long", 1);
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  myPL.set<float>("my float", 4.3);
  myPL.set<double>("my double", 4.3);
  myPL.set("my string", "hello");
  myPL.set("my char", 'c');
  myPL.set("my bool", true);

  // Array are also supported for the following types
  myPL.set<Array<int> >("my array int", tuple<int>(1, 2));
  myPL.set<Array<unsigned int> >("my array unsigned int",
    tuple<unsigned int>(1));
  myPL.set<Array<short int> > ("my array short int",
    tuple<short int>(1, 2));
  myPL.set<Array<unsigned short int> > ("my array unsigned short int",
    tuple<unsigned short int>(1, 2));
  myPL.set<Array<long int> >("my array long int",
    tuple<long int>(1, 2));
  myPL.set<Array<unsigned long int> >("my array unsigned long int",
    tuple<unsigned long int>(1, 2));
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  myPL.set<Array<long long int> >("my array long long int",
    tuple<long long int>(1, 2));
  myPL.set<Array<unsigned long long int> >("my array unsigned long long int",
    tuple<unsigned long long>(1, 2));
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  myPL.set<Array<float> >("my array float", tuple<float>(1,1, 2.2));
  myPL.set<Array<double> >("my array double", tuple<double>(1,1, 2.2));
  myPL.set<Array<std::string> >("my array string",
    tuple<std::string>("hello", "world"));

  //Now for the custom data type. First, lets put one in the parameter list.
  CustomDataType sampleCustom(3, "hello");

  myPL.set("my custom data", sampleCustom);
  
  //Now before we write this out to xml, we have to make sure we have a
  //converter for our cusomt data type. Since our custom datatype overrides
  //the operator<< and operator>> we can just use and instance of the
  //StandardTemplatedParameterConverter.  We'll do this using the convience
  //macro. Look at the source code for the macro to see everything that's
  //actually goiing on. It's in Teuchos_ParameterEntryXMLConverterDB.hpp.

  TEUCHOS_ADD_TYPE_CONVERTER(CustomDataType);

  //Now we'll write it out to xml.
  Teuchos::writeParameterListToXmlFile(myPL, "xml_data_types_test_list.xml");
  //Then read it in to a new list.

  Teuchos::writeParameterListToXmlOStream(
    myPL,
    std::cout);
 
  const RCP<ParameterList> readPL = 
    Teuchos::getParametersFromXmlFile("xml_data_types_test_list.xml");

  std::cout << *readPL;

  //If we compare them, we'll see they're equal
  if(*readPL == myPL){
    std::cout << "Huzzah!\n";
  }
  else{
    std::cerr << "Uh oh...";
    return -1;
  }

  // Read the parameters in one at a time
  const int myInt = readPL->get<int>("my int");
  std::cout << "myInt = " << myInt << "\n";
  const unsigned int myUnsignedInt = readPL->get<unsigned int>("my unsigned int");
  std::cout << "myUnsignedInt = " << myUnsignedInt << "\n";
  const short int myShortInt = readPL->get<short int>("my short int");
  std::cout << "myShortInt = " << myShortInt << "\n";
  const short int myShort = readPL->get<short>("my short");
  std::cout << "myShort = " << myShort << "\n";
  const unsigned short int myUnsignedShortInt = readPL->get<unsigned short int>("my unsigned short int");
  std::cout << "myUnsignedShortInt = " << myUnsignedShortInt << "\n";
  const unsigned short int myUnsignedShort = readPL->get<unsigned short>("my unsigned short");
  std::cout << "myUnsignedShort = " << myUnsignedShort << "\n";
  const long int myLongInt = readPL->get<long int>("my long int");
  std::cout << "myLongInt = " << myLongInt << "\n";
  const long int myLong = readPL->get<long>("my long");
  std::cout << "myLong = " << myLong << "\n";
  const unsigned long int myUnsignedLongInt = readPL->get<unsigned long int>("my unsigned long int");
  std::cout << "myUnsignedLongInt = " << myUnsignedLongInt << "\n";
  const unsigned long myUnsignedLong = readPL->get<unsigned long>("my unsigned long");
  std::cout << "myUnsignedLong = " << myUnsignedLong << "\n";
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  const long long int myLongLongInt = readPL->get<long long int>("my long long int");
  std::cout << "myLongLongInt = " << myLongLongInt << "\n";
  const long long int myLongLong = readPL->get<long long>("my long long");
  std::cout << "myLongLong = " << myLongLong << "\n";
  const unsigned long long int myUnsignedLongLongInt = readPL->get<unsigned long long int>("my unsigned long long int");
  std::cout << "myUnsignedLongLongInt = " << myUnsignedLongLongInt << "\n";
  const unsigned long long myUnsignedLongLong = readPL->get<unsigned long long>("my unsigned long long");
  std::cout << "myUnsignedLongLong = " << myUnsignedLongLong << "\n";
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  const float myFloat = readPL->get<float>("my float");
  std::cout << "myFloat = " << myFloat << "\n";
  const double myDouble = readPL->get<double>("my double");
  std::cout << "myDouble = " << myDouble << "\n";
  const std::string myString = readPL->get<std::string>("my string");
  std::cout << "myString = " << myString << "\n";
  const char myChar = readPL->get<char>("my char");
  std::cout << "myChar = " << myChar << "\n";
  const bool myBool = readPL->get<bool>("my bool");
  std::cout << "myBool = " << myBool << "\n";
  const Array<int> myIntArray = readPL->get<Array<int> >("my array int");
  std::cout << "myIntArray = " << myIntArray << "\n";
  const Array<float> myFloatArray = readPL->get<Array<float> >("my array float");
  std::cout << "myFloatArray = " << myFloatArray << "\n";
  const Array<double> myDoubleArray = readPL->get<Array<double> >("my array double");
  std::cout << "myDoubleArray = " << myDoubleArray << "\n";
  const Array<std::string> myStringArray = readPL->get<Array<std::string> >("my array string");
  std::cout << "myStringArray = " << myStringArray << "\n";

  /**
   * Final Notes: StandardTemplatedParameterConverter should suit most your
   * needs. Buf if for some reason you don't feel like overrideing the
   * inseration and extraction operators, you can allways subclass the
   * ParameterEntryXMLConverter class and do your own thing.
   */
  return 0;
}
