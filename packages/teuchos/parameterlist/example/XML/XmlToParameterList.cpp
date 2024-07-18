// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Version.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

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
