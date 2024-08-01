// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_epetra_stratimikos_solver.hpp"
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
      = Thyra::test_epetra_stratimikos_solver(
        &paramList, dumpAll, verbose?&*out:0
        );

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success)

  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
