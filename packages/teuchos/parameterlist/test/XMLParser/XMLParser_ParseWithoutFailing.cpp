// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_UnitTestHarness.hpp"

std::string filename;
std::string error;

namespace Teuchos {

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::UnitTestRepository::getCLP().setOption(
      "filename", &filename, "XML file to parse" );
  }

  TEUCHOS_UNIT_TEST( XMLParser, IgnoreDeclarations )
  {
    FileInputSource src(filename);
    Teuchos::RCP< Teuchos::XMLInputStream > stream = src.stream();
    XMLParser parser(stream);
    bool caughtError = false;
    try {
      Teuchos::XMLObject obj = parser.parse();
      out << "Parsed XML object: \n" << obj.toString() << std::endl;
    }
    catch (std::runtime_error &err) {
      caughtError = true;
      out << "XML Parser caught exception:\n" << err.what() << "\n";
    }
    TEST_EQUALITY_CONST(caughtError, false);
  }

} // namespace Teuchos
