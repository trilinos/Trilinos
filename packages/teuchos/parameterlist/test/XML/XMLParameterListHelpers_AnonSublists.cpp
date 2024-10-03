// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
