// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


namespace {


// Helper function to compute a single norm for a vector
double epetraNorm2( const Epetra_Vector &v )
{
  double norm[1] = { -1.0 };
  v.Norm2(&norm[0]);
  return norm[0];
}


} // namespace


int main(int argc, char* argv[])
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::CommandLineProcessor;
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;

  bool success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {


    //
    // A) Program setup code
    //

    //
    // Read options from command-line
    //

    std::string     matrixFile             = "";
    std::string     extraParamsFile        = "";
    double          tol                    = 1e-5;
    bool            onlyPrintOptions       = false;
    bool            printXmlFormat         = false;
    bool            showDoc                = true;

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    CommandLineProcessor  clp(false); // Don't throw exceptions

    // Set up command-line options for the linear solver that will be used!
    linearSolverBuilder.setupCLP(&clp);

    clp.setOption( "matrix-file", &matrixFile
                   ,"Defines the matrix and perhaps the RHS and LHS for a linear system to be solved." );
    clp.setOption( "tol", &tol, "Tolerance to check against the scaled residual norm." );
    clp.setOption( "extra-params-file", &extraParamsFile, "File containing extra parameters in XML format.");
    clp.setOption( "only-print-options", "continue-after-printing-options", &onlyPrintOptions
                   ,"Only print options and stop or continue on" );
    clp.setOption( "print-xml-format", "print-readable-format", &printXmlFormat
                   ,"Print the valid options in XML format or in readable format." );
    clp.setOption( "show-doc", "hide-doc", &showDoc
                   ,"Print the valid options with or without documentation." );

    clp.setDocString(
      "Simple example for the use of the Stratimikos facade Stratimikos::DefaultLinearSolverBuilder.\n"
      "\n"
      "To print out just the valid options use --matrix-file=\"\" --only-print-options with --print-xml-format"
      " or --print-readable-format.\n"
      "\n"
      "To solve a linear system read from a file use --matrix-file=\"SomeFile.mtx\""
      " with options read from an XML file using --linear-solver-params-file=\"SomeFile.xml\"\n"
      );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Print out the valid options and stop if asked
    //

    if(onlyPrintOptions) {
      if(printXmlFormat)
        Teuchos::writeParameterListToXmlOStream(
          *linearSolverBuilder.getValidParameters()
          ,*out
          );
      else
        linearSolverBuilder.getValidParameters()->print(
          *out,PLPrintOptions().indent(2).showTypes(true).showDoc(showDoc)
          );
      return 0;
    }


    //
    // B) Epetra-specific code that sets up the linear system to be solved
    //
    // While the below code reads in the Epetra objects from a file, you can
    // setup the Epetra objects any way you would like.  Note that this next
    // set of code as nothing to do with Thyra at all, and it should not.
    //

    *out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\' ...\n";

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    RCP<Epetra_CrsMatrix> epetra_A;
    RCP<Epetra_Vector> epetra_x, epetra_b;
    EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A, NULL, &epetra_x, &epetra_b );

    if(!epetra_b.get()) {
      *out << "\nThe RHS b was not read in so generate a new random vector ...\n";
      epetra_b = rcp(new Epetra_Vector(epetra_A->OperatorRangeMap()));
      epetra_b->Random();
    }

    if(!epetra_x.get()) {
      *out << "\nThe LHS x was not read in so generate a new zero vector ...\n";
      epetra_x = rcp(new Epetra_Vector(epetra_A->OperatorDomainMap()));
      epetra_x->PutScalar(0.0); // Initial guess is critical!
    }

    *out << "\nPrinting statistics of the Epetra linear system ...\n";

    *out
      << "\n  Epetra_CrsMatrix epetra_A of dimension "
      << epetra_A->OperatorRangeMap().NumGlobalElements()
      << " x " << epetra_A->OperatorDomainMap().NumGlobalElements()
      << "\n  ||epetraA||inf = " << epetra_A->NormInf()
      << "\n  ||epetra_b||2 = " << epetraNorm2(*epetra_b)
      << "\n  ||epetra_x||2 = " << epetraNorm2(*epetra_x)
      << "\n";


    //
    // C) The "Glue" code that takes Epetra objects and wraps them as Thyra
    // objects
    //
    // This next set of code wraps the Epetra objects that define the linear
    // system to be solved as Thyra objects so that they can be passed to the
    // linear solver.
    //

    RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp( epetra_A );
    RCP<Thyra::VectorBase<double> > x = Thyra::create_Vector( epetra_x, A->domain() );
    RCP<const Thyra::VectorBase<double> > b = Thyra::create_Vector( epetra_b, A->range() );

    // Note that above Thyra is only interacted with in the most trival of
    // ways.  For most users, Thyra will only be seen as a thin wrapper that
    // they need know little about in order to wrap their objects in order to
    // pass them to Thyra-enabled solvers.


    //
    // D) Thyra-specific code for solving the linear system
    //
    // Note that this code has no mention of any concrete implementation and
    // therefore can be used in any use case.
    //

    // Reading in the solver parameters from the parameters file and/or from
    // the command line.  This was setup by the command-line options
    // set by the setupCLP(...) function above.
    linearSolverBuilder.readParameters(out.get());

    // Augment parameters if appropriate
    if(extraParamsFile.length()) {
      Teuchos::updateParametersFromXmlFile( "./"+extraParamsFile,
        linearSolverBuilder.getNonconstParameterList().ptr() );
    }

    // Create a linear solver factory given information read from the
    // parameter list.
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
      linearSolverBuilder.createLinearSolveStrategy("");

    // Setup output stream and the verbosity level
    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

    // Create a linear solver based on the forward operator A
    RCP<Thyra::LinearOpWithSolveBase<double> > lows =
      Thyra::linearOpWithSolve(*lowsFactory, A);

    // Solve the linear system (note: the initial guess in 'x' is critical)
    Thyra::SolveStatus<double> status =
      Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
    *out << "\nSolve status:\n" << status;

    // Write the linear solver parameters after they were read
    linearSolverBuilder.writeParamsFile(*lowsFactory);


    //
    // E) Post process the solution and check the error
    //
    // Note that the below code is based only on the Epetra objects themselves
    // and does not in any way depend or interact with any Thyra-based
    // objects.  The point is that most users of Thyra can largely gloss over
    // the fact that Thyra is really being used for anything.
    //

    // Wipe out the Thyra wrapper for x to guarantee that the solution will be
    // written back to epetra_x!  At the time of this writting this is not
    // really needed but the behavior may change at some point so this is a
    // good idea.
    x = Teuchos::null;

    *out
      << "\nSolution ||epetra_x||2 = " << epetraNorm2(*epetra_x) << "\n";

    *out << "\nTesting the solution error ||b-A*x||/||b|| computed through the Epetra objects ...\n";

    // r = b - A*x
    Epetra_Vector epetra_r(*epetra_b);
    {
      Epetra_Vector epetra_A_x(epetra_A->OperatorRangeMap());
      epetra_A->Apply(*epetra_x,epetra_A_x);
      epetra_r.Update(-1.0,epetra_A_x,1.0);
    }

    const double
      nrm_r = epetraNorm2(epetra_r),
      nrm_b = epetraNorm2(*epetra_b),
      rel_err = ( nrm_r / nrm_b );
    const bool
      passed = (rel_err <= tol);

    *out
      << "||b-A*x||/||b|| = " << nrm_r << "/" << nrm_b << " = " << rel_err
      << " < tol = " << tol << " ? " << ( passed ? "passed" : "failed" ) << "\n";

    if(!passed) success = false;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)

  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

}
