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

#include <Teuchos_Parser.hpp>

using Teuchos::ParameterList;
using Teuchos::RCP;

namespace TeuchosTests
{
  TEUCHOS_UNIT_TEST(YAML, XmlEquivalence)
  {
    using std::string;
    using std::vector;
    vector<string> matchStems;
    matchStems.push_back("Match1");
    matchStems.push_back("Match2");
    matchStems.push_back("Match3");
    matchStems.push_back("Match4");
    for(size_t i = 0; i < matchStems.size(); i++)
    {
      string yamlFile = matchStems[i] + ".yaml";
      RCP<ParameterList> yamlList = Teuchos::getParametersFromYamlFile(yamlFile);
      string xmlFile =  matchStems[i] + ".xml";
      RCP<ParameterList> xmlList = Teuchos::getParametersFromXmlFile(xmlFile);
      TEST_EQUALITY(Teuchos::haveSameValues(*xmlList, *yamlList, true), true);
    }
  }
  TEUCHOS_UNIT_TEST(YAML, IntVsDouble)
  {
    //YAML2 has a double param that has the same name/value as an int param in XML2
    //YAML reader should recognize the double and the param lists should not be equivalent
    RCP<ParameterList> xmlList = Teuchos::getParametersFromXmlFile("IntVsDouble.xml");
    RCP<ParameterList> yamlList = Teuchos::getParametersFromYamlFile("IntVsDouble.yaml");
    TEST_EQUALITY(Teuchos::haveSameValues(*xmlList, *yamlList), false);
  }
  TEUCHOS_UNIT_TEST(YAML, IllegalKeyString)
  {
    TEST_THROW(Teuchos::getParametersFromYamlFile("IllegalKeyString.yaml");, Teuchos::ParserFail);
  }

  TEUCHOS_UNIT_TEST(YAML, Issue1801)
  {
    Teuchos::getParametersFromYamlString(
        "My Awesome Problem:\n"
        "  Particle Periodic:\n"
        "    X: \"-1.0, 1.0\"\n"
        "  emotions: happy_sad, indifferent\n"
        "...\n"
        );
  }

  TEUCHOS_UNIT_TEST(YAML, PR1805)
  {
    RCP<ParameterList> params = Teuchos::getParametersFromYamlString(
        "My Awesome Problem:\n"
        "\tMesh:\n"
        "\t\tInline:\n"
        "\t\t\tType: Quad\n"
        "\t\t\tElements: [     10,     10 ]\n"
        "...\n"
        );
  }

  TEUCHOS_UNIT_TEST(YAML, IntAndDoubleArray)
  {
    int correctInts[5] = {2, 3, 5, 7, 11};
    //the last number with 10 dec digits of precision should test correct double conversion
    double correctDoubles[5] = {2.718, 3.14159, 1.618, 1.23456789, 42.1337};
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
  }
  TEUCHOS_UNIT_TEST(YAML, InconsistentArrayType)
  {
    std::string correctStrings[5] = {"2", "3", "5", "7", "imastring"};
    double correctDoubles[5] = {2, 3, 1.618, 1.23456789, 42.1337};
    RCP<ParameterList> plist = Teuchos::getParametersFromYamlFile("InconsistentArray.yaml");
    //verify that stringArray and doubleArray have the correct types and the correct values
    const Teuchos::Array<std::string>& stringArr = plist->get<Teuchos::Array<std::string> >("stringArray");
    const Teuchos::Array<double>& doubleArr = plist->get<Teuchos::Array<double> >("doubleArray");
    for(int i = 0; i < 5; i++)
    {
      if(stringArr[i] != correctStrings[i])
      {
        throw std::runtime_error(std::string("stringArray value is incorrect."));
      }
      if(doubleArr[i] != correctDoubles[i])
      {
        throw std::runtime_error(std::string("doubleArray value is incorrect."));
      }
    }
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
      "      <Parameter name=\"Bool-looking String\" type=\"string\" value=\"TRUE\"/>\n"
      "    </ParameterList>\n"
      "  </ParameterList>\n";
    RCP<ParameterList> xmlParams = Teuchos::getParametersFromXmlString(xmlString);
    std::stringstream yamlOutStream;
    yamlOutStream << std::showpoint << std::fixed << std::setprecision(1);
    Teuchos::YAMLParameterList::writeYamlStream(yamlOutStream, *xmlParams);
    std::string yamlString = yamlOutStream.str();
    std::string expectedYamlString =
      "%YAML 1.1\n"
      "---\n"
      "ANONYMOUS:\n"
      "  Problem: \n"
      "    Neumann BCs: \n"
      "      Time Dependent NBC on SS cyl_outside for DOF all set P: \n"
      "        BC Values: [[0.0], [10.0], [20.0]]\n"
      "  Discretization: \n"
      "    Node Set Associations: [['1', '2'], [top, bottom]]\n"
      "    Bool-looking String: 'TRUE'\n"
      "...\n";
    TEST_EQUALITY(yamlString, expectedYamlString);
    std::stringstream yamlInStream(yamlString);
    RCP<ParameterList> yamlParams;
    yamlParams = Teuchos::YAMLParameterList::parseYamlStream(yamlInStream);
    std::stringstream yamlOutStream2;
    yamlOutStream2 << std::showpoint << std::fixed << std::setprecision(1);
    Teuchos::YAMLParameterList::writeYamlStream(yamlOutStream2, *yamlParams);
    std::string yamlString2 = yamlOutStream2.str();
    TEST_EQUALITY(yamlString2, expectedYamlString);
  }

  TEUCHOS_UNIT_TEST(YAML, Issue1815)
  {
    Teuchos::getParametersFromYamlString(
      "Header:\n"
      "  Output:\n"
      "    File Name: electrostatic.exo\n"
      "    Cell Average Quantities:\n"
      "      eblock-0_0: ES_POTENTIAL, E0\n"
      "      \n"
      "  Particle Dump:\n"
      "    File Name: beam_emit2.h5part\n");
  }

  TEUCHOS_UNIT_TEST(YAML, Issue1807part2)
  {
    Teuchos::getParametersFromYamlString(
      "Header:\n"
      "  Particle Dump:\n"
      "    File Name: beam_emit2.h5part\n"
      "#   Stride Time: 5.0e-12\n");
  }

  TEUCHOS_UNIT_TEST(YAML, Issue2090)
  {
    auto params = Teuchos::getParametersFromYamlString(
      "Parameter List:\n"
      "  Boundary Conditions:\n"
      "    Bottom:\n"
      "      Dirichlet:\n"
      "        Sideset: bottom\n"
      "        Field:   ES_POTENTIAL\n"
      "        Value: |\n"
      "          double r_sq = xin*xin+yin*yin;\n"
      "          double factor = 0.5*1.e8*1.60217662e-19/(2*3.14159265358979323846*8.854187817e-12);\n"
      "          ES_POTENTIAL= factor*log(r_sq) +3*xin-3*yin;\n"
      "  # end Boundary Conditions\n");
    TEST_EQUALITY(
        Teuchos::getParameter<std::string>(
           params->sublist("Boundary Conditions", true)
           .sublist("Bottom", true)
           .sublist("Dirichlet", true)
          ,"Value"),
      "double r_sq = xin*xin+yin*yin;\n"
      "double factor = 0.5*1.e8*1.60217662e-19/(2*3.14159265358979323846*8.854187817e-12);\n"
      "ES_POTENTIAL= factor*log(r_sq) +3*xin-3*yin;\n");
  }

  TEUCHOS_UNIT_TEST(YAML, Issue2306)
  {
    // ensure that duplicate names throw an exception
    TEST_THROW(Teuchos::getParametersFromYamlString(
      "Foo:\n"
      "  Bar:\n"
      "    Value: 1\n"
      "  Bar:\n"
      "    Value: 2\n"),
      Teuchos::ParserFail);
  }

} //namespace TeuchosTests

