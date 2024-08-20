// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_single_belos_thyra_solver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


int main(int argc, char* argv[])
{
  
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //
    
    std::string     matrixFile             = "";
    bool            testTranspose          = false;
    bool            usePreconditioner      = true;
    int             numRhs                 = 1;
    int             numRandomVectors       = 1;
    double          maxFwdError            = 1e-14;
    int             maxIterations          = 400;
    int             maxRestarts            = 25;
    int             gmresKrylovLength      = 25;
    int             outputFrequency        = 10;
    bool            outputMaxResOnly       = true;
    int             blockSize              = 1;
    double          maxResid               = 1e-6;
    double          maxSolutionError       = 1e-6;
    bool            showAllTests           = false;
    bool            dumpAll                = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "matrix-file", &matrixFile, "Matrix input file [Required]." );
    clp.setOption( "test-transpose", "no-test-transpose", &testTranspose, "Test the transpose solve or not." );
    clp.setOption( "use-preconditioner", "no-use-preconditioner", &usePreconditioner, "Use the preconditioner or not." );
    clp.setOption( "num-rhs", &numRhs, "Number of RHS in linear solve." );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "max-fwd-error", &maxFwdError, "The maximum relative error in the forward operator." );
    clp.setOption( "max-iters", &maxIterations, "The maximum number of linear solver iterations to take." );
    clp.setOption( "max-restarts", &maxRestarts, "???." );
    clp.setOption( "gmres-krylov-length", &gmresKrylovLength, "???." );
    clp.setOption( "output-frequency", &outputFrequency, "Number of linear solver iterations between output" );
    clp.setOption( "output-max-res-only", "output-all-res", &outputMaxResOnly, "Determines if only the max residual is printed or if all residuals are printed per iteration." );
    clp.setOption( "block-size", &blockSize, "???." );
    clp.setOption( "max-resid", &maxResid, "The maximum relative error in the residual." );
    clp.setOption( "max-solution-error", &maxSolutionError, "The maximum relative error in the solution of the linear system." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEUCHOS_TEST_FOR_EXCEPT( matrixFile == "" );

    Teuchos::ParameterList belosLOWSFPL;

    belosLOWSFPL.set("Solver Type","Block GMRES");

    Teuchos::ParameterList& belosLOWSFPL_solver =
      belosLOWSFPL.sublist("Solver Types");

    Teuchos::ParameterList& belosLOWSFPL_gmres =
      belosLOWSFPL_solver.sublist("Block GMRES");

    belosLOWSFPL_gmres.set("Maximum Iterations",int(maxIterations));
    belosLOWSFPL_gmres.set("Convergence Tolerance",double(maxResid));
    belosLOWSFPL_gmres.set("Maximum Restarts",int(maxRestarts));
    belosLOWSFPL_gmres.set("Block Size",int(blockSize));
    belosLOWSFPL_gmres.set("Num Blocks",int(gmresKrylovLength));
    belosLOWSFPL_gmres.set("Output Frequency",int(outputFrequency));
    belosLOWSFPL_gmres.set("Show Maximum Residual Norm Only",bool(outputMaxResOnly));

    Teuchos::ParameterList precPL("Ifpack");
    if(usePreconditioner) {
      precPL.set("Overlap",int(2));
      precPL.set("Prec Type","ILUT");
    }
    
    success
      = Thyra::test_single_belos_thyra_solver(
        matrixFile,testTranspose,usePreconditioner,numRhs,numRandomVectors
        ,maxFwdError,maxResid,maxSolutionError,showAllTests,dumpAll
        ,&belosLOWSFPL,&precPL
        ,verbose?&*out:0
        );

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? 0 : 1 );
}
