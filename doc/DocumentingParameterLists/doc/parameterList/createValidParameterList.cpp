// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_VerboseObject.hpp"
#include <fstream>

// Add #include statements here to pull in the code that will be
// needed to generate one or more validated ParameterLists

// Flatten the namespaces of certain classes
using std::string;
using Teuchos::ParameterList;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;

int main(int argc, char* argv[])
{

  bool success = false;
  bool verbose = false;
  try
  {
    //
    // Open an output file stream for writing our XML file
    //
    std::ofstream out;
    out.open("stratimikos.xml", std::ofstream::out);

    //
    // We will print to the 'out' ofstream, and that output will be
    // valid XML describing the validated ParameterList.  For the
    // purposes of generating nicely-formatted HTML documentation for
    // this ParameterList, we also need to include an XSL header line.
    // This bool will control whether we include this header line,
    // which can be controlled at the command line.
    //
    bool xsl_header_flag = true;

    //
    // Set up the command line processor.  All versions of this
    // executable should support the add-xsl-header /
    // suppress-xsl-header command line options.  If you want a single
    // executable to support multiple ParameterLists, you could put
    // additional options here to control which ParameterList to
    // output.
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
    // Here is where code should go that is required to generate the
    // validated XML file.  If your class uses a construct-then-init
    // idiom, then this is where the default constructor would be
    // called.  If this executable supports ParameterLists for more
    // than one class, then this is where the logic to pick the
    // appropriate class would go.
    //

    //
    // Write the validated ParameterList to the fancy output stream.
    // Replace PARAMETERLIST with a reference to the
    // Teuchos::ParameterList obtained in the step above.
    //
    Teuchos::writeParameterListToXmlOStream(PARAMETERLIST, out);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

}
