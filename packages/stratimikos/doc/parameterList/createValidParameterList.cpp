// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include <fstream>

// Flatten the namespaces of certain classes
using std::string;
using Teuchos::ParameterList;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;

int main(int argc, char* argv[])
{

  bool success = false;
  bool verbose = false;
  try {
    std::ofstream out;
    out.open("stratimikos.xml", std::ofstream::out);

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    //
    // We will print to standard out, and that output will be valid XML
    // describing the validated ParameterList.  For the purposes of
    // generating nicely-formatted HTML documentation for this
    // ParameterList, we also need to include an XSL header line.  This
    // bool will control whether we include this header line, which can
    // be controlled at the command line.
    //
    bool xsl_header_flag = true;

    //
    // Set up the command line processor.  All versions of this
    // executable should support the add-xsl-header/suppress-xsl-header
    // command line options.  If you want a single executable to support
    // multiple ParameterLists, you could put additional options here to
    // control which ParameterList to output.
    //
    CommandLineProcessor clp(false);  //don't throw exceptions
    clp.recogniseAllOptions(true);
    clp.setOption("add-xsl-header",
      "suppress-xsl-header",
      &xsl_header_flag, 
      "XSL header flag");

    //
    // Parse the command line and quit if not successful
    //
    CommandLineProcessor::EParseCommandLineReturn parse_return =
      clp.parse(argc, argv);
    if(parse_return != CommandLineProcessor::PARSE_SUCCESSFUL)
      return parse_return;

    //
    // Print the XSL header line if requested
    //
    if (xsl_header_flag )
      out << "<?xml-stylesheet type=\"text/xsl\" "
          << "href=\"common/parameterList/parameterList.xsl\"?>\n";

    //
    // Obtain the validated ParameterList and write it to the fancy
    // output stream.  If you wanted to support multiple ParameterLists,
    // this is where the logic would go to choose between them.  Note
    // that Domi has a function that returns the valid ParameterList,
    // but that a more common use case will be to construct a class
    // (that supports the construct-then-init paradigm) with a default
    // constructor and then call its getValidParameters() method.
    //

    Teuchos::writeParameterListToXmlOStream(*linearSolverBuilder.getValidParameters(), out);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
