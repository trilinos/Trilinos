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

using std::string;
using std::endl;

namespace Teuchos {


/* Simple test of Teuchos XMLParser class */
TEUCHOS_UNIT_TEST( XMLParser, orderedWriteRead )
{
  out << endl;

  /* create a ParameterList object */
  string xmlstring1;
  ParameterList wrotepl1("Parent");
  {
    ParameterList c1("Child1");
    ParameterList c2("Child2");
    c1.set("cp1", "first1");
    c1.set("cp2", "second1");
    c2.set("cp3", "first2");
    c2.set("cp4", "second2");
    wrotepl1.set("FirstSublist",c1);
    wrotepl1.set("SecondSublist",c2);
    /* create an XML object from the ParameterList and write it to a string */
    XMLParameterListWriter xml2pl;
    XMLObject xmlprob = xml2pl.toXML(wrotepl1);
    std::ostringstream ss;
    ss << xmlprob;
    xmlstring1 = ss.str();
    out << "*** String 1" << endl;
    out << xmlstring1 << endl;
  }

  string xmlstring2;
  ParameterList wrotepl2("Parent");
  {
    ParameterList c1("Child1");
    ParameterList c2("Child2");
    // swap the ordering
    c1.set("cp2", "second1");
    c1.set("cp1", "first1");
    c2.set("cp4", "second2");
    c2.set("cp3", "first2");
    wrotepl2.set("SecondSublist",c2);
    wrotepl2.set("FirstSublist",c1);
    /* create an XML object from the ParameterList and write it to a string */
    XMLParameterListWriter xml2pl;
    XMLObject xmlprob = xml2pl.toXML(wrotepl2);
    std::ostringstream ss;
    ss << xmlprob;
    xmlstring2 = ss.str();
    out << "*** String 2" << endl;
    out << xmlstring2 << endl;
  }

  // the different PL orderings should be reflected in the ParameterLists and their XML string representations
  TEST_INEQUALITY(wrotepl1, wrotepl2);
  TEST_INEQUALITY(xmlstring1, xmlstring2);

  /* create a input source, parser to read the string */
  ParameterList readpl1, readpl2;
  {
    StringInputSource src(xmlstring1);
    XMLParser parser(src.stream());
    XMLObject xmlprob = parser.parse();
    XMLParameterListReader pl2xml;
    readpl1 = pl2xml.toParameterList(xmlprob);
  }
  {
    StringInputSource src(xmlstring2);
    XMLParser parser(src.stream());
    XMLObject xmlprob = parser.parse();
    XMLParameterListReader pl2xml;
    readpl2 = pl2xml.toParameterList(xmlprob);
  }

  /* check that the parameter lists do not match */
  TEST_INEQUALITY(readpl1, readpl2);

}


TEUCHOS_UNIT_TEST( XMLParser, simpleOrderedRead )
{

  /* create a ParameterList object */
  string xmlstring1;
  ParameterList plGold("ParentList");
  {
    ParameterList c1("Z");
    ParameterList c2("A");
    c1.set("A", "first1");
    c1.set("Z", "second1");
    c2.set("Z", "first2");
    c2.set("A", "second2");
    plGold.set("9FirstSublist",c1);
    plGold.set("1SecondSublist",c2);
  }

  string xmlsrc(
    "<ParameterList name=\"ParentList\">\n"
    "  <ParameterList name=\"9FirstSublist\">\n"
    "    <Parameter name=\"A\" type=\"string\" value=\"first1\"/>\n"
    "    <Parameter name=\"Z\" type=\"string\" value=\"second1\"/>\n"
    "  </ParameterList>\n"
    "  <ParameterList name=\"1SecondSublist\">\n"
    "    <Parameter name=\"Z\" type=\"string\" value=\"first2\"/>\n"
    "    <Parameter name=\"A\" type=\"string\" value=\"second2\"/>\n"
    "  </ParameterList>\n"
    "</ParameterList>\n");
  ParameterList plTest;
  {
    StringInputSource src(xmlsrc);
    XMLParser parser(src.stream());
    XMLObject xmlprob = parser.parse();
    XMLParameterListReader pl2xml;
    plTest = pl2xml.toParameterList(xmlprob);
  }

  TEST_EQUALITY(plTest, plGold);
}


} // namespace Teuchos
