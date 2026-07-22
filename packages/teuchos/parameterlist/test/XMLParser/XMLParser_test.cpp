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

namespace Teuchos {


/* Simple test of Teuchos XMLParser class */
TEUCHOS_UNIT_TEST( XMLParser, simpletest )
{

  /* create a ParameterList object */
  ParameterList problem("Problem");
  ParameterList solver("Solver");
  ParameterList prec("Preconditioner");
  prec.set("type", "ILUk");           // set some of these to isUsed for completeness
  TEST_NOTHROW( prec.get<std::string>("type") );
  prec.set("k", 2);
  solver.set("Preconditioner",prec);
  solver.set("type", "gmres");        // set some of these to isUsed for completeness
  TEST_NOTHROW( solver.get<std::string>("type") );
  solver.set("maxiters", 1000);
  solver.set("restarts", 100);
  solver.set("special1","\"&\"");     // test the XML outputting and parsing for correctness
  solver.set("special2","\'&\'");     // test the XML outputting and parsing for correctness
  solver.set("special3","\"&\'");     // test the XML outputting and parsing for correctness
  solver.set("tol", 1.0e-10);         // set some of these to isUsed for completeness
  TEST_NOTHROW( solver.get<double>("tol") );
  problem.set("Solver",solver);

  std::cout << "*** ParameterList (original)" << std::endl;
  std::cout << problem << std::endl;

  /* create an XML object from the ParameterList */
  XMLParameterListWriter xml2pl;
  XMLObject xmlprob1 = xml2pl.toXML(problem);

  /* write the XML to a std::string */
  std::ostringstream ss;
  ss << xmlprob1;
  std::string strproblem = ss.str();

  std::cout << "*** XML from ParameterListParameterListWriter.toXML().toString()" << std::endl;
  std::cout << xmlprob1 << std::endl;

  /* create a input source, parser to read the std::string */
  StringInputSource src(strproblem);
  XMLParser parser(src.stream());

  /* parse XML in a std::string */
  XMLObject xmlprob2 = parser.parse();

  std::cout << "*** XML from XMLParser.parse()" << std::endl;
  std::cout << xmlprob2 << std::endl;

  /* convert the new XML object to a ParameterList */
  XMLParameterListReader pl2xml;
  ParameterList problem2 = pl2xml.toParameterList(xmlprob2);

  std::cout << "*** ParameterList from XMLParameterListReader.toParameterList()" << std::endl;
  std::cout << problem2 << std::endl;

  /* check that the new parameter list matches the old one */
  TEST_EQUALITY(problem2, problem);

}


} // namespace Teuchos



