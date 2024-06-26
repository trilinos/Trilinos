// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************// @HEADER

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
#include "Tpetra_MultiVector.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"

typedef Tpetra::MultiVector<> multivector_type;
typedef Tpetra::Operator<> operator_type;
typedef multivector_type::scalar_type scalar_type;
typedef Belos::SolverManager<scalar_type, multivector_type, operator_type> solver_type;

// Flatten the namespaces of certain classes
using std::string;
using Teuchos::ParameterList;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::parameterList;

void writeXSLHeader( std::ofstream& out )
{
      out << "<?xml-stylesheet type=\"text/xsl\" "
          << "href=\"common/parameterList/parameterList.xsl\"?>\n";
}

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
    // Here is where code should go that is required to generate the
    // validated XML file.  If your class uses a construct-then-init
    // idiom, then this is where the default constructor would be
    // called.  If this executable supports ParameterLists for more
    // than one class, then this is where the logic to pick the
    // appropriate class would go.
    //

	Belos::SolverFactory<scalar_type, multivector_type, operator_type> sFactory;
	RCP<solver_type> solver;
	for ( int i = 0; i < 4; ++i )
	{
		RCP<ParameterList> nullParams = parameterList();
		std::cout << "writing parameter list " << i+1 << " of 9." << std::endl;
		if ( i == 0 )
		{
			solver = sFactory.create("Block GMRES", nullParams);
			out.open("belos_BlockGmres.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> gmresParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *gmresParams, out);
			out.close();
		}
		else if ( i == 1 )
		{
			solver = sFactory.create("Pseudo Block GMRES", nullParams);
			out.open("belos_PseudoBlockGmres.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> pseudoGmresParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *pseudoGmresParams, out);
			out.close();
		}
		else if ( i == 2 )
		{
			solver = sFactory.create("Block CG", nullParams);
			out.open("belos_BlockCG.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> blockCgParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *blockCgParams, out);
			out.close();
		}
		else if ( i == 3 )
		{
			solver = sFactory.create("Pseudo Block CG", nullParams);
			out.open("belos_PseudoBlockCG.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> pseudoBlockCgParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *pseudoBlockCgParams, out);
			out.close();
		}
	}

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

}
