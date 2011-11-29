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

#include "Teuchos_Version.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <fstream>

int main( int argc, char* argv[] )
{
  
  using Teuchos::inoutArg;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  std::cout << std::endl << Teuchos::Teuchos_Version() << std::endl;

  bool success = true;

  try {

    std::string    xmlInFileName = "";
    std::string    extraXmlFile = "";
    std::string    xmlOutFileName = "paramList.out";

    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setOption("extra-xml-file",&extraXmlFile,"File with extra XML text that will modify the initial XML read in");
    clp.setOption("xml-out-file",&xmlOutFileName,"The XML file to write the final parameter list to");
    clp.setDocString(
      "This example program shows how to read in a parameter list from an"
      " XML file (given by --xml-in-file=xmlInFileName) and then modify it"
      " given some XML specified on the command-line (given by --extra-xml=extrXmlStr)."
      " The final parameter list is then written back to an XML file."
      " (given by --xml-out-file=xmlOutFileName)."
      );
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
      std::cout << "\nEnd Result: TEST FAILED" << std::endl;
      return parse_return;
    }

    Teuchos::ParameterList paramList;

    if(xmlInFileName.length()) {
      std::cout << "\nReading a parameter list from the XML file \""<<xmlInFileName<<"\" ...\n";
      Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(paramList));
      std::cout << "\nParameter list read from the XML file \""<<xmlInFileName<<"\":\n\n";
      paramList.print(std::cout,2,true,true);
    }
    
    std::string line("");
    if(extraXmlFile.length()) {
      std::ifstream myfile(extraXmlFile.c_str());
      if (myfile.is_open())
      {
        getline (myfile,line);
        std::cout << line << "\n";
        myfile.close();
      }
      std::cout << "\nUpdating the parameter list given the extra XML std::string:\n\n"<<line<<"\n";
      Teuchos::updateParametersFromXmlString(line, inoutArg(paramList));
      std::cout << "\nParameter list after ammending extra XML std::string:\n\n";
      paramList.print(std::cout,2,true,true);
    }

    std::cout << "\nWriting the final parameter list back to the XML file \""<<xmlOutFileName<<"\" ... \n";
    Teuchos::writeParameterListToXmlFile(paramList,xmlOutFileName);

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

  if(success)
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  else
    std::cout << "\nEnd Result: TEST FAILED" << std::endl;

  return ( success ? 0 : 1 );

}
