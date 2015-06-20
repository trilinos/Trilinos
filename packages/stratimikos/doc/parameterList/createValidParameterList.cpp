// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
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
