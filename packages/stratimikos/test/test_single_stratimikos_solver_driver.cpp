// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "test_single_stratimikos_solver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //
    
    std::string     inputFile              = "";
    std::string     extraParams            = "";
    bool            dumpAll                = false;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
    clp.setOption( "input-file", &inputFile, "Input file [Required].", true );
    clp.setOption( "extra-params", &extraParams, "Extra parameters overriding the parameters read in from --input-file");
    clp.setDocString(
      "Testing program for Trilinos (and non-Trilinos) linear solvers access through Thyra."
      );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    Teuchos::ParameterList paramList;
    if(verbose) *out << "\nReading parameters from XML file \""<<inputFile<<"\" ...\n";
    Teuchos::updateParametersFromXmlFile(inputFile, Teuchos::inOutArg(paramList));
    if(extraParams.length()) {
      if(verbose) *out << "\nAppending extra parameters from the XML string \""<<extraParams<<"\" ...\n";
      Teuchos::updateParametersFromXmlString(extraParams, Teuchos::inOutArg(paramList));
    }
    
    success
      = Thyra::test_single_stratimikos_solver(
        &paramList, dumpAll, verbose?&*out:0
        );
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success)
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? 0 : 1 );
}
