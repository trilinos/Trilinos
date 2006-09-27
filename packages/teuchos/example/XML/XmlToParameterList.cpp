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

#include "Teuchos_Version.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  std::cout << std::endl << Teuchos::Teuchos_Version() << std::endl;

  bool success = true;

  try {

    std::string    xmlInFileName = "";
    std::string    extraXmlStr = "";
    std::string    xmlOutFileName = "paramList.out";

    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setOption("extra-xml",&extraXmlStr,"Extra XML that will modify the initial XML read in");
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
      cout << "\nEnd Result: TEST FAILED" << endl;
      return parse_return;
    }

    Teuchos::ParameterList paramList;

    if(xmlInFileName.length()) {
      std::cout << "\nReading a parameter list from the XML file \""<<xmlInFileName<<"\" ...\n";
      Teuchos::updateParametersFromXmlFile(xmlInFileName,&paramList);
      std::cout << "\nParameter list read from the XML file \""<<xmlInFileName<<"\":\n\n";
      paramList.print(std::cout,2,true,true);
    }
    
    if(extraXmlStr.length()) {
      std::cout << "\nUpdating the parameter list given the extra XML string:\n\n"<<extraXmlStr<<"\n";
      Teuchos::updateParametersFromXmlString(extraXmlStr,&paramList);
      std::cout << "\nParameter list after ammending extra XML string:\n\n";
      paramList.print(std::cout,2,true,true);
    }

    std::cout << "\nWriting the final parameter list back to the XML file \""<<xmlOutFileName<<"\" ... \n";
    Teuchos::writeParameterListToXmlFile(paramList,xmlOutFileName);

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

  if(success)
    std::cout << "\nEnd Result: TEST PASSED" << endl;
  else
    std::cout << "\nEnd Result: TEST FAILED" << endl;

  return ( success ? 0 : 1 );

}
