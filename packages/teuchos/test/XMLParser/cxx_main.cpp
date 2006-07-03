// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListReader.hpp"

using namespace Teuchos;
using std::string;

/* Test of Teuchos XMLParser class */

int main(int argc, void** argv)
{
  cout << Teuchos::Teuchos_Version() << endl << endl;

  bool testfailed = false;

  try
    {
      MPISession::init(&argc, &argv);

      /* create a ParameterList object */
      ParameterList problem("Problem");
      ParameterList solver("Solver");
      ParameterList prec("Preconditioner");
      prec.set("type", "ILUk");
      prec.set("k", 2);
      solver.set("Preconditioner",prec);
      solver.set("type", "gmres");
      solver.set("maxiters", 1000);
      solver.set("restarts", 100);
      solver.set("tol", 1.0e-10);
      problem.set("Solver",solver);

      cout << "*** ParameterList" << endl;
      cout << problem << endl;

      /* create an XML object from the ParameterList */
      XMLParameterListWriter xml2pl;
      XMLObject xmlprob1 = xml2pl.toXML(problem);

      /* write the XML to a string */
      string strproblem = xmlprob1.toString();

      cout << "*** ParameterList.toXML" << endl;
      cout << xmlprob1 << endl;

      /* create a input source, parser to read the string */
      StringInputSource src(strproblem);
      XMLParser parser(src.stream());

      /* parse XML in a string */
      XMLObject xmlprob2 = parser.parse();

      cout << "*** parser.parse" << endl;
      cout << xmlprob2 << endl;

      /* convert the new XML object to a ParameterList */
      XMLParameterListReader pl2xml;
      ParameterList problem2 = pl2xml.toParameterList(xmlprob2);

      cout << "*** XML.toParameterList" << endl;
      cout << problem2 << endl;

      /* check that the new parameter list matches the old one */
      if (problem2 != problem) {
        testfailed = true;
      }
    }
  catch(std::exception& e)
    {
      cerr << e.what() << endl;
      testfailed = true;
    }
  MPISession::finalize();

  if (testfailed) {
    cout << "End Result: TEST FAILED" << endl;
    return -1;
  }

  cout << "End Result: TEST PASSED" << endl;
  return 0;
}
