// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_YamlParser_decl.hpp>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <fstream>
#include <iomanip>
#include <sstream>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::DefaultComm;

namespace TeuchosTests
{
  TEUCHOS_UNIT_TEST(YAML, MPIBroadcast)
  {
    //load Match1.xml and Match1.yaml on proc 0, broadcast them, and make sure it matches on all procs
    //note: the following code is run on all procs
    RCP<ParameterList> xmlList = rcp(new ParameterList);
    RCP<ParameterList> yamlList = rcp(new ParameterList);

    RCP<const Teuchos::Comm<int> > comm = DefaultComm<int>::getComm ();
    const std::string n1("Match1.xml");
    Teuchos::updateParametersFromXmlFileAndBroadcast(n1, xmlList.ptr(), *comm, true);
    const std::string n2("Match1.yaml");
    Teuchos::updateParametersFromYamlFileAndBroadcast(n2, yamlList.ptr(), *comm, true);
    TEST_EQUALITY(Teuchos::haveSameValues(*xmlList, *yamlList), true);
  }
  TEUCHOS_UNIT_TEST(YAML, ConvertFromXML)
  {
    using std::string;

    RCP<const Teuchos::Comm<int> > comm = DefaultComm<int>::getComm ();

    //This list can contain any valid XML param lists in the unit_tests/yaml/
    std::vector<string> xmlFiles;
    xmlFiles.push_back("Match1.xml");
    xmlFiles.push_back("Match2.xml");
    xmlFiles.push_back("Match3.xml");
    xmlFiles.push_back("Match4.xml");
    xmlFiles.push_back("input_restingHydrostatic_RK4.xml");
    xmlFiles.push_back("plasma_oscillation_rtc.xml");
    for(size_t i = 0; i < xmlFiles.size(); i++)
    {
      //reading from XML
      std::ifstream xmlStream(xmlFiles[i].c_str());
      //emitting converted YAML
      std::ostringstream yamlStream;
      yamlStream << std::setprecision(17) << std::scientific;
      Teuchos::convertXmlToYaml(xmlStream, yamlStream);
      //now read back both formats to compare
      RCP<ParameterList> xmlList = Teuchos::getParametersFromXmlFile(xmlFiles[i]);
      string yamlText = yamlStream.str();
      string debugYamlFileName = xmlFiles[i] + ".yaml";
      {
      std::ofstream debugYamlFileStream(debugYamlFileName.c_str());
      debugYamlFileStream << yamlText;
      }
      RCP<ParameterList> yamlList = Teuchos::YAMLParameterList::parseYamlText(yamlText,
          debugYamlFileName);
      TEST_EQUALITY(Teuchos::haveSameValues(*xmlList, *yamlList, true), true);
    }
  }

} //namespace TeuchosTests

