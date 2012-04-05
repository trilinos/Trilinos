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
