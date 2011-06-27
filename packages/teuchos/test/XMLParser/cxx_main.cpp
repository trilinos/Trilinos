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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListReader.hpp"

using std::string;
using Teuchos::ParameterList;
using Teuchos::XMLObject;
using Teuchos::XMLParser;
using Teuchos::XMLParameterListReader;
using Teuchos::XMLParameterListWriter;
using Teuchos::StringInputSource;

/* Test of Teuchos XMLParser class */

int main(int argc, char** argv)
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  bool testfailed = false;
  double tmp;
      
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try
    {

      /* create a ParameterList object */
      ParameterList problem("Problem");
      ParameterList solver("Solver");
      ParameterList prec("Preconditioner");
      prec.set("type", "ILUk");           // set some of these to isUsed for completeness
      {
        std::string l_tmp = prec.get<std::string>("type");
      }
      prec.set("k", 2);
      solver.set("Preconditioner",prec);
      solver.set("type", "gmres");        // set some of these to isUsed for completeness
      {
        std::string l_tmp = solver.get<std::string>("type");
      }
      solver.set("maxiters", 1000);
      solver.set("restarts", 100);
      solver.set("special1","\"&\"");     // test the XML outputting and parsing for correctness
      solver.set("special2","\'&\'");     // test the XML outputting and parsing for correctness
      solver.set("special3","\"&\'");     // test the XML outputting and parsing for correctness
      solver.set("tol", 1.0e-10);         // set some of these to isUsed for completeness
      {
        tmp = solver.get<double>("tol");
      }
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
      if (problem2 != problem) {
        testfailed = true;
      }
    }
  catch(std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      testfailed = true;
    }

  if (testfailed) {
    std::cout << "End Result: TEST FAILED" << std::endl;
    return -1;
  }

  std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}
