// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  TEUCHOS_UNIT_TEST(YAML, PR1805) // "\t" has been removed since yaml does not support tabs
  {
    RCP<ParameterList> params = Teuchos::getParametersFromYamlString(
        "My Awesome Problem:\n"
        "  Mesh:\n"
        "    Inline:\n"
        "      Type: Quad\n"
        "      Elements: [     10,     10 ]\n"
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
      "      <Parameter name=\"Bool-looking String\" type=\"string\" value=\"TRUE\" docString=\"my docString\"/>\n"
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

  TEUCHOS_UNIT_TEST(YAML, keep_top_name)
  {
    Teuchos::ParameterList pl;
    char const * const cstr =
      "%YAML 1.1\n"
      "---\n"
      "Albany:\n"
      "  some param: 5\n"
      "...\n";
    Teuchos::updateParametersFromYamlCString(cstr, Teuchos::ptr(&pl), true);
    std::stringstream ss;
    ss << std::showpoint;
    Teuchos::writeParameterListToYamlOStream(pl, ss);
    auto s = ss.str();
    TEST_EQUALITY(s, cstr);
  }

  TEUCHOS_UNIT_TEST(YAML, long_long_param)
  {
    auto pl = Teuchos::getParametersFromYamlString(
      "List:\n"
      " small number: 54\n"
      " small double: 3.0\n"
      " scientific: 3.123e02\n"
      " big number: 72057594037927936\n"
    );
    TEST_EQUALITY(pl->isType<int>("small number"), true);
    TEST_EQUALITY(pl->isType<double>("small double"), true);
    TEST_EQUALITY(pl->isType<double>("scientific"), true);
    TEST_EQUALITY(pl->isType<long long>("big number"), true);
    TEST_EQUALITY(pl->get<long long>("big number"), 72057594037927936ll);
  }

  TEUCHOS_UNIT_TEST(YAML, bools)
  {
    auto pl = Teuchos::getParametersFromYamlString(
      "List:\n"
      " input_true: true\n"
      " input_false: false\n"
      " input_TRUE: TRUE\n"
      " input_FALSE: FALSE\n"
      " input_True: True\n"
      " input_False: False\n"
      " input_yes: yes\n"
      " input_no: no\n"
    );
    TEST_EQUALITY(pl->isType<bool>("input_true"), true);
    TEST_EQUALITY(pl->isType<bool>("input_false"), true);
    TEST_EQUALITY(pl->isType<bool>("input_yes"), true);
    TEST_EQUALITY(pl->isType<bool>("input_no"), true);
    TEST_EQUALITY(pl->isType<bool>("input_TRUE"), true);
    TEST_EQUALITY(pl->isType<bool>("input_True"), true);
    TEST_EQUALITY(pl->isType<bool>("input_FALSE"), true);
    TEST_EQUALITY(pl->isType<bool>("input_False"), true);
    TEST_EQUALITY(pl->get<bool>("input_true"), true);
    TEST_EQUALITY(pl->get<bool>("input_false"), false);
    TEST_EQUALITY(pl->get<bool>("input_yes"), true);
    TEST_EQUALITY(pl->get<bool>("input_no"), false);
    TEST_EQUALITY(pl->get<bool>("input_TRUE"), true);
    TEST_EQUALITY(pl->get<bool>("input_True"), true);
    TEST_EQUALITY(pl->get<bool>("input_FALSE"), false);
    TEST_EQUALITY(pl->get<bool>("input_False"), false);
  }

  TEUCHOS_UNIT_TEST(YAML, flow_map)
  {
    auto pl = Teuchos::getParametersFromYamlString(
      "List:\n"
      " Fields: {rho: 0.125, px: 0., py: 0., pz: 0., rho_E: 0.25}\n");
    auto& field_pl = pl->sublist("Fields");
    TEST_EQUALITY(field_pl.get<double>("rho"), 0.125);
  }

  TEUCHOS_UNIT_TEST(YAML, root_name)
  {
    Teuchos::ParameterList pl;
    Teuchos::updateParametersFromYamlString(
      "mycode:\n"
      "  sublist:\n"
      "    param1: foo\n",
      Teuchos::ptr(&pl),
      true,
      "root_name test"
      );
    auto& sublist = pl.sublist("sublist");
    TEST_EQUALITY(sublist.name(), "mycode->sublist");
  }

  TEUCHOS_UNIT_TEST(YAML, null_node)
  {
    RCP<ParameterList> pl = Teuchos::getParametersFromYamlString(
        "mycode:\n"
        "  empty_node:\n"
    );
    TEST_EQUALITY(pl->isSublist("empty_node"), true);
  }

#ifdef HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
  TEUCHOS_UNIT_TEST(YAML, yamlcpp_parser)
  {
    RCP<ParameterList> pl = Teuchos::getParametersFromYamlString(
        "mycode:\n"
        "  list_of_2d_arrays:\n"
        "    - [[1,2,3], [4,5,6]]\n"
        "    - [[7,8,9], [10,11,12]]\n"
        "  ragged_array:\n"
        "    - [1,2,3]\n"
        "    - [1,2,3,4]\n"
        "  line_continuation: [\n"
        "    1,2,3,\n"
        "    4,5,6\n"
        "  ]\n"
        "  # allow unicode comments: ±\n"
    );

    using threeDarr_t = Teuchos::Array<Teuchos::Array<Teuchos::Array<int>>>;
    threeDarr_t& list_of_arrs = pl->get<threeDarr_t>("list_of_2d_arrays");
    threeDarr_t correct_list_of_arrs = {
      {{1, 2, 3}, {4, 5, 6}},
      {{7, 8, 9}, {10, 11, 12}}
    };
    for (int i=0; i<list_of_arrs.size(); i++) {
      for (int j=0; j<list_of_arrs[i].size(); j++) {
        for (int k=0; k<list_of_arrs[i][j].size(); k++) {
          TEST_EQUALITY(correct_list_of_arrs[i][j][k], list_of_arrs[i][j][k]);
        }
      }
    }

    using twoDarr_t = Teuchos::Array<Teuchos::Array<int>>;
    twoDarr_t ragged_arr = pl->get<twoDarr_t>("ragged_array");
    twoDarr_t correct_ragged_arr = {
      {1, 2, 3},
      {1, 2, 3, 4}
    };
    for (int i=0; i<ragged_arr.size(); i++) {
      for (int j=0; j<ragged_arr[i].size(); j++) {
        TEST_EQUALITY(correct_ragged_arr[i][j], ragged_arr[i][j]);
      }
    }

    using arr_t   = Teuchos::Array<int>;
    arr_t arr = pl->get<arr_t>("line_continuation");
    arr_t correct_arr = {1, 2, 3, 4, 5, 6};
    for (int i=0; i<arr.size(); i++) {
      TEST_EQUALITY(correct_arr[i], arr[i]);
    }
  }

  TEUCHOS_UNIT_TEST(YAML, yaml_throws)
  {
  TEST_THROW(Teuchos::getParametersFromYamlString(
    "Foo:\n"
    "  [60,2,3]: 1\n"),
    Teuchos::YamlKeyError)
  TEST_THROW(Teuchos::getParametersFromYamlString(
    "Foo:\n"
    "  Array:\n"
    "  - 1.3e0.2\n"
    "  - [1,2,3]"),
    Teuchos::YamlSequenceError)
  TEST_THROW(Teuchos::getParametersFromYamlString(
    "Foo: 1\n"),
    Teuchos::YamlStructureError)
  }
  // It is not clear how to test Teuchos::YamlUndefinedNode, but the throw
  // is left in the source code to protect against any unforeseen cases.

#endif // HAVE_TEUCHOSPARAMETERLIST_YAMLCPP

} //namespace TeuchosTests
