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

#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"

using std::string;
using std::endl;

namespace Teuchos {

  TEUCHOS_UNIT_TEST( XMLParameterListHelpers, anonymousSublists )
  {
    // the sublists are missing names; that is allowed now
    string anon_sublist(
      "<ParameterList name=\"Parent\">\n"
      "  <ParameterList>\n"
      "    <Parameter name=\"param_a\" type=\"string\" value=\"a\"/>\n"
      "    <ParameterList>\n"
      "      <Parameter name=\"param_b\" type=\"string\" value=\"b\"/>\n"
      "    </ParameterList>\n"
      "  </ParameterList>\n"
      "</ParameterList>\n");
    StringInputSource src(anon_sublist);
    XMLParser parser(src.stream());
    XMLObject xmlprob = parser.parse();
    XMLParameterListReader pl2xml;
    ParameterList plTest;
    TEST_NOTHROW( plTest = pl2xml.toParameterList(xmlprob) );
    out << plTest;
  }

  TEUCHOS_UNIT_TEST( XMLParameterListHelpers, anonymousParam )
  {
    // one parameter is missing a name; but all parameters are required to have a name
    string anon_param(
      "<ParameterList name=\"Parent\">\n"
      "  <ParameterList>\n"
      "    <Parameter name=\"param_a\" type=\"string\" value=\"a\"/>\n"
      "    <ParameterList>\n"
      "      <Parameter type=\"string\" value=\"b\"/>\n"
      "    </ParameterList>\n"
      "  </ParameterList>\n"
      "</ParameterList>\n");
    StringInputSource src(anon_param);
    XMLParser parser(src.stream());
    XMLObject xmlprob = parser.parse();
    XMLParameterListReader pl2xml;
    ParameterList plTest;
    // this exception name is mis-spelled; oh well
    TEST_THROW( plTest = pl2xml.toParameterList(xmlprob), NoNameAttributeExecption);
  }

} // namespace Teuchos
