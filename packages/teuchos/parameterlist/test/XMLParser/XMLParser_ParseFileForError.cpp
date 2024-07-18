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

  TEUCHOS_UNIT_TEST( XMLParser, ParseFile )
  {
    FileInputSource src(filename);
    Teuchos::RCP< Teuchos::XMLInputStream > stream = src.stream();
    while (1) {
      unsigned char c;
      TEUCHOS_TEST_FOR_EXCEPTION( stream->readBytes(&c,1) < 1, std::runtime_error, "Failure reading error message from test file." );
      if (c == '\n') break;
      if (c != '\r') error.push_back(c);
    }
    out << "Expected error string: \"" << error << "\"" << std::endl;
    bool caughtError = false;
    try {
      XMLParser parser(stream);
      parser.parse();
    }
    catch (std::runtime_error &err) {
      std::string what = err.what();
      caughtError = true;
      out << "Caught exception:\n" << what << "\n";
      TEST_INEQUALITY_CONST( what.find(error), std::string::npos );
    }
    TEST_EQUALITY_CONST(caughtError, true);
  }

} // namespace Teuchos
