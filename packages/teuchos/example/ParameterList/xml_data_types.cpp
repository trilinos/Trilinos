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
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  Teuchos::ParameterList My_List;

  //Basic data types
  My_List.set("int", 1);
  My_List.set("unsigned int", (unsigned int)1);
  My_List.set("short int", (short)1);
  My_List.set("unsigned short int", (unsigned short)1);
  My_List.set("long int", (long)1);
  My_List.set("unsigned long int", (unsigned long)1);
  #ifdef HAVE_TEUCHOS_LONG_LONG_INT
  My_List.set("long long int", (long long)1);
  My_List.set("unsigned long long int", (unsigned long long)1);
  #endif
  My_List.set("float", (float)4.3);
  My_List.set("double", (double)4.3);
  My_List.set("string", "hello");
  My_List.set("char", 'c');
  My_List.set("bool", true);

  //Teuchos::Array are also supported for the following types
  Teuchos::Array<int> intArray = Teuchos::tuple<int>(1);
  My_List.set("array int", intArray );
  Teuchos::Array<unsigned int> uintArray = Teuchos::tuple<unsigned int>(1);
  My_List.set("array unsigned int", uintArray);
  Teuchos::Array<short> shortArray = Teuchos::tuple<short>(1);
  My_List.set("array short int", shortArray);
  Teuchos::Array<unsigned short> ushortArray = Teuchos::tuple<unsigned short>(1);
  My_List.set("array unsigned short int", ushortArray);
  Teuchos::Array<long> longArray = Teuchos::tuple<long>(1);
  My_List.set("array long int", longArray);
  Teuchos::Array<unsigned long> ulongArray = Teuchos::tuple<unsigned long>(1);
  My_List.set("array unsigned long int", ulongArray);
  #ifdef HAVE_TEUCHOS_LONG_LONG_INT
  Teuchos::Array<long long> longlongArray = Teuchos::tuple<long long>(1);
  My_List.set("array long long int", longlongArray);
  Teuchos::Array<unsigned long long> ulonglongArray = Teuchos::tuple<unsigned long long>(1);
  My_List.set("array unsigned long long int", ulonglongArray);
  #endif
  Teuchos::Array<float> floatArray = Teuchos::tuple<float>(4.3);
  My_List.set("array float", floatArray);
  Teuchos::Array<double> doubleArray = Teuchos::tuple<double>(4.3);
  My_List.set("array double", doubleArray);
  Teuchos::Array<std::string> stringArray = Teuchos::tuple<std::string>("hello");
  My_List.set("array string", stringArray);

  //Now for the custom data type. First, lets put one in the parameter list.
  CustomDataType sampleCustom(3, "hello");

  My_List.set("custom data", sampleCustom);
  
  //Now before we write this out to xml, we have to make sure we have a converter 
  //for our cusomt data type. Since our custom datatype overrides the operator<< and
  //operator>> we can just use and instance of the StandardTemplatedParameterConverter.
  //We'll do this using the convience macro. Look at the source code for the macro
  //to see everything that's actually goiing on. It's in Teuchos_ParameterEntryXMLConverterDB.hpp .

  TEUCHOS_ADD_TYPE_CONVERTER(CustomDataType);

  //Now we'll write it out to xml.
  Teuchos::writeParameterListToXmlFile(My_List, "xml_data_types_test_list.xml");
  //Then read it in to a new list.

  Teuchos::writeParameterListToXmlOStream(
    My_List,
    std::cout);

 
  Teuchos::RCP<Teuchos::ParameterList> readIn = Teuchos::getParametersFromXmlFile("xml_data_types_test_list.xml");

  std::cout << *readIn;
  //If we compare them, we'll see they're equal
  if(*readIn == My_List){
    std::cout << "Huzzah!\n";
  }
  else{
    std::cerr << "Uh oh...";
    return -1;
  }

  /**
   * Final Notes:
   * StandardTemplatedParameterConverter should suit most your needs. Buf if for some reason you
   * don't feel like overrideing the inseration and extraction operators, you can allways subclass
   * the ParameterEntryXMLConverter class and do your own thing.
   */
  return 0;
}
