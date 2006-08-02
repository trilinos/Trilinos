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

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char* argv[])
{
  using Teuchos::rcp;
  using Teuchos::RefCountPtr;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;
  using Teuchos::getParameter;
  using Teuchos::CommandLineProcessor;

  bool result, success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //
    
    std::string     inputFile              = "";
    std::string     extraParams            = "";
    bool            onlyPrintOptions       = false;
    bool            printXmlFormat         = false;

    Thyra::DefaultRealLinearSolverBuilder linearSolverBuilder;

    CommandLineProcessor  clp(false); // Don't throw exceptions

    linearSolverBuilder.setupCLP(&clp);

    clp.setOption( "only-print-options", "continue-after-printing-options", &onlyPrintOptions
                   ,"Only print options and stop or continue on" );
    clp.setOption( "print-xml-format", "print-readable-format", &printXmlFormat
                   ,"Print the valid options in XML format or in readable format." );
                   
    clp.setDocString(
      "Simple example for the use of Stratimikos software.  Right now this only prints"
      " out valid options and then stops.\n\n"
      "To print out just the valid options use --input-file=\"\" --only-print-options with --print-xml-format"
      " or --print-readable-format."
      );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Print out the valid options if asked to
    //

    if(onlyPrintOptions) {
      if(printXmlFormat)
        Teuchos::writeParameterListToXmlOStream(
          *linearSolverBuilder.getValidParameters()
          ,*out
          );
      else
        linearSolverBuilder.getValidParameters()->print(*out,1,true,false);
      return 0;
    }

    //
    // Reading in the solver parameters
    //

    linearSolverBuilder.readParameters(out.get());

    //
    // Solve a linear system or something ...
    //

    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

    // ToDo: Solve a linear system or something!

    //
    // Write the linear solver parameter after they where read
    //
    
    linearSolverBuilder.writeParamsFile(*lowsFactory);
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success)
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? 0 : 1 );
}
