// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_YamlParser_decl.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include <fstream>

using namespace Teuchos;

int main(int argc, char *argv[]) {

  TEUCHOS_ASSERT(argc == 2);

  std::string file_prefix(argv[1]);

  std::string input_xml_file_name(file_prefix+".xml");
  auto pList = getParametersFromXmlFile(input_xml_file_name);

  std::string output_yaml_file_name(file_prefix+".yaml");
  writeParameterListToYamlFile(*pList,output_yaml_file_name);

  return 0;
}
