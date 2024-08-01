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
