// @HEADER
//
// ***********************************************************************
//
//         Teuchos
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_YamlParser_decl.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_YamlParameterListCoreHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_TwoDArray.hpp>



using Teuchos::ParameterList;
using Teuchos::RCP;

namespace TeuchosTests
{
  TEUCHOS_UNIT_TEST(YAML, XmlEquivalence)
  {
    using std::string;
    using std::vector;
    vector<string> matchStems = {"Match1", "Match2", "Match3", "Match4"};
    for(size_t i = 0; i < matchStems.size(); i++)
    {
      string xmlFile =  matchStems[i] + ".xml";
      string yamlFile = matchStems[i] + ".yaml";
      RCP<ParameterList> xmlList = Teuchos::getParametersFromXmlFile(xmlFile);
      RCP<ParameterList> yamlList = Teuchos::getParametersFromYamlFile(yamlFile);
      TEST_EQUALITY(Teuchos::haveSameValuesUnordered(*xmlList, *yamlList), true);
    }
  }
  TEUCHOS_UNIT_TEST(YAML, IntVsDouble)
  {
    //YAML2 has a double param that has the same name/value as an int param in XML2
    //YAML reader should recognize the double and the param lists should not be equivalent
    RCP<ParameterList> xmlList = Teuchos::getParametersFromXmlFile("IntVsDouble.xml");
    RCP<ParameterList> yamlList = Teuchos::getParametersFromYamlFile("IntVsDouble.yaml");
    TEST_EQUALITY(Teuchos::haveSameValuesUnordered(*xmlList, *yamlList), false);
  }
  TEUCHOS_UNIT_TEST(YAML, IllegalKeyString)
  {
    TEST_THROW(Teuchos::getParametersFromYamlFile("IllegalKeyString.yaml");, YAML::ParserException);
  }
  TEUCHOS_UNIT_TEST(YAML, IntAndDoubleArray)
  {
    int correctInts[5] = {2, 3, 5, 7, 11};
    //the last number with 10 dec digits of precision should test correct double conversion
    double correctDoubles[5] = {2.718, 3.14159, 1.618, 1.23456789, 42.1337};
    TEST_NOTHROW({
      RCP<Teuchos::ParameterList> params = Teuchos::getParametersFromYamlFile("Arrays.yaml");
      //Retrieve arrays from a specific sublist (tests the mixed nesting of sequence/map)
      ParameterList& sublist = params->get<ParameterList>("smoother: params");
      Teuchos::Array<int>& intArr = sublist.get<Teuchos::Array<int> >("intArray");
      Teuchos::Array<double>& doubleArr = sublist.get<Teuchos::Array<double> >("doubleArray");
      TEST_EQUALITY(intArr.size(), 5);
      TEST_EQUALITY(doubleArr.size(), 5);
      for(int i = 0; i < 5; i++)
      {
        TEST_EQUALITY(correctInts[i], intArr[i]);
        TEST_EQUALITY(correctDoubles[i], doubleArr[i]);
      }
    });
  }
  TEUCHOS_UNIT_TEST(YAML, InconsistentArrayType)
  {
    std::string correctStrings[5] = {"2", "3", "5", "7", "imastring"};
    double correctDoubles[5] = {2, 3, 1.618, 1.23456789, 42.1337};
    TEST_NOTHROW({
      RCP<ParameterList> plist = Teuchos::getParametersFromYamlFile("InconsistentArray.yaml");
      //verify that stringArray and doubleArray have the correct types and the correct values
      const Teuchos::Array<std::string>& stringArr = plist->get<Teuchos::Array<std::string> >("stringArray");
      const Teuchos::Array<double>& doubleArr = plist->get<Teuchos::Array<double> >("doubleArray");
      for(int i = 0; i < 5; i++)
      {
        if(stringArr[i] != correctStrings[i])
        {
          throw std::runtime_error(std::string("stringArray[") + std::to_string(i) + "] is incorrect.");
        }
        if(doubleArr[i] != correctDoubles[i])
        {
          throw std::runtime_error(std::string("doubleArray value [") + std::to_string(i) + "] is incorrect.");
        }
      }
    });
  }
  TEUCHOS_UNIT_TEST(YAML, TwoDArrayConvert)
  {
    std::string xmlString =
      "  <ParameterList>\n"
      "    <ParameterList name=\"Problem\">\n"
      "      <ParameterList name=\"Neumann BCs\">\n"
      "        <ParameterList name=\"Time Dependent NBC on SS cyl_outside for DOF all set P\">\n"
      "          <Parameter name=\"BC Values\" type=\"TwoDArray(double)\" value=\"3x1:{ 0.0, 10.0, 20.0}\"/>\n"
      "        </ParameterList>\n"
      "      </ParameterList>\n"
      "    </ParameterList>\n"
      "    <ParameterList name=\"Discretization\">\n"
      "      <Parameter name=\"Node Set Associations\" type=\"TwoDArray(string)\" value=\"2x2:{1, 2, top, bottom}\"/>\n"
      "    </ParameterList>\n"
      "  </ParameterList>\n";
    RCP<ParameterList> xmlParams = Teuchos::getParametersFromXmlString(xmlString);
    std::stringstream yamlOutStream;
    Teuchos::YAMLParameterList::writeYamlStream(yamlOutStream, *xmlParams);
    std::string yamlString = yamlOutStream.str();
    std::string expectedYamlString =
      "%YAML 1.1\n"
      "---\n"
      "ANONYMOUS:\n"
      "  Problem: \n"
      "    Neumann BCs: \n"
      "      Time Dependent NBC on SS cyl_outside for DOF all set P: \n"
      "        BC Values: [[0], [10], [20]]\n"
      "  Discretization: \n"
      "    Node Set Associations: [['1', '2'], [top, bottom]]\n"
      "...\n";
    TEST_EQUALITY(yamlString, expectedYamlString);
    std::stringstream yamlInStream(yamlString);
    RCP<ParameterList> yamlParams;
    yamlParams = Teuchos::YAMLParameterList::parseYamlStream(yamlInStream);
    std::stringstream yamlOutStream2;
    Teuchos::YAMLParameterList::writeYamlStream(yamlOutStream2, *yamlParams);
    std::string yamlString2 = yamlOutStream2.str();
  /* There are issues with older versions of yaml-cpp not maintaining the order
     of parameters in a list. see Trilinos issue #1268.
     that is why we use the Unordered comparison instead of the pure text one. */
  //TEST_EQUALITY(yamlString2, expectedYamlString);
    std::stringstream yamlInStream2(yamlString2);
    RCP<ParameterList> yamlParams2;
    yamlParams2 = Teuchos::YAMLParameterList::parseYamlStream(yamlInStream2);
    TEST_EQUALITY(Teuchos::haveSameValuesUnordered(*yamlParams, *yamlParams2), true);
  }
} //namespace TeuchosTests

