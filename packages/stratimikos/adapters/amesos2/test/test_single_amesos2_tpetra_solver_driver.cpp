// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_single_amesos2_tpetra_solver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Tpetra_Core.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


int main(int argc, char* argv[])
{

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  auto comm = Tpetra::getDefaultComm();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //

    std::string     matrixFile             = "";
    std::string     solverType             = "";
    int             numRhs                 = 1;
    int             numRandomVectors       = 1;
    double          maxFwdError            = 1e-14;
    int             outputFrequency        = 10;
    bool            outputMaxResOnly       = true;
    double          maxResid               = 1e-6;
    double          maxSolutionError       = 1e-6;
    bool            showAllTests           = false;
    bool            dumpAll                = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "matrix-file", &matrixFile, "Matrix input file [Required]." );
    clp.setOption( "solver-type", &solverType, "Type of solver to use [Required]." );
    clp.setOption( "num-rhs", &numRhs, "Number of RHS in linear solve." );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "max-fwd-error", &maxFwdError, "The maximum relative error in the forward operator." );
    clp.setOption( "output-frequency", &outputFrequency, "Number of linear solver iterations between output" );
    clp.setOption( "output-max-res-only", "output-all-res", &outputMaxResOnly, "Determines if only the max residual is printed or if all residuals are printed per iteration." );
    clp.setOption( "max-resid", &maxResid, "The maximum relative error in the residual." );
    clp.setOption( "max-solution-error", &maxSolutionError, "The maximum relative error in the solution of the linear system." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEUCHOS_TEST_FOR_EXCEPT( matrixFile == "" );
    TEUCHOS_TEST_FOR_EXCEPT( solverType == "" );

    Teuchos::ParameterList amesos2LOWSFPL("Amesos2");

    amesos2LOWSFPL.set("Solver Type",solverType);

    success
      = Thyra::test_single_amesos2_tpetra_solver(
        matrixFile,numRhs,numRandomVectors
        ,maxFwdError,maxResid,maxSolutionError,showAllTests,dumpAll
        ,&amesos2LOWSFPL
        ,verbose?&*out:0
        ,comm
        );

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)

  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? 0 : 1 );
}
