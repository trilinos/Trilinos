// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_single_aztecoo_thyra_solver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

int main(int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::FancyOStream out(Teuchos::rcp(&std::cout,false));

  try {

    //
    // Read options from command-line
    //
    
    std::string     matrixFile             = "";
    bool            testTranspose          = true;
    int             numRandomVectors       = 1;
    double          maxFwdError            = 1e-14;
    int             maxIterations          = 400;
    double          maxResid               = 1e-6;
    double          maxSolutionError       = 1e-6;
    bool            showAllTests           = false;
    bool            dumpAll                = false;
    std::string     aztecOutputLevel       = "freq";
    int             aztecOutputFreq        = 0;
    std::string     aztecPrec              = "none";
    std::string     aztecSubdomainSolve    = "ilu";

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "matrix-file", &matrixFile, "Matrix input file [Required]." );
    clp.setOption( "test-transpose", "no-test-transpose", &testTranspose, "Test the transpose solve or not." );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "max-fwd-error", &maxFwdError, "The maximum relative error in the forward operator." );
    clp.setOption( "max-iters", &maxIterations, "The maximum number of linear solver iterations to take." );
    clp.setOption( "max-resid", &maxResid, "The maximum relative error in the residual." );
    clp.setOption( "max-solution-error", &maxSolutionError, "The maximum relative error in the solution of the linear system." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    clp.setOption( "aztec-output-level", &aztecOutputLevel, "Aztec output level (freq,last,summary,warnings,all)" );
    clp.setOption( "aztec-output-freq", &aztecOutputFreq, "Aztec output freqency (> 0)" );
    clp.setOption( "aztec-prec", &aztecPrec, "Type of aztec preconditioner (none,sym_GS,Neumann,Jacobi,ls,dom_decomp)" );
    clp.setOption( "aztec-subdomain-solve", &aztecSubdomainSolve, "Type of subdomain solve for --aztec-prec==dom_decomp only (ilu,ilut)" );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEUCHOS_TEST_FOR_EXCEPT( matrixFile == "" );

    Teuchos::ParameterList
      paramList("AztecOOLinearOpWithSolve");
    Teuchos::ParameterList
      &fwdSolveParamList = paramList.sublist("Forward Solve"),
      &adjSolveParamList = paramList.sublist("Adjoint Solve");
    fwdSolveParamList.set("Max Iterations",maxIterations);
    adjSolveParamList.set("Max Iterations",maxIterations);
/*
    Teuchos::ParameterList
      &fwdAztecOOParamList = fwdSolveParamList.sublist("AztecOO"),
      &adjAztecOOParamList = fwdSolveParamList.sublist("AztecOO");
    if( aztecOutputLevel != "freq" ) {
      fwdAztecOOParamList.set("Output Frequency",aztecOutputLevel);
      adjAztecOOParamList.set("Output Frequency",aztecOutputLevel);
    }
    else {
      fwdAztecOOParamList.set("Output Frequency",aztecOutputFreq);
      adjAztecOOParamList.set("Output Frequency",aztecOutputFreq);
    }
    if( aztecPrec != "none" ) {
      fwdAztecOOParamList.set("Aztec Preconditioner",aztecPrec);
      adjAztecOOParamList.set("Aztec Preconditioner",aztecPrec);
      if(aztecPrec=="dom_decomp") {
        fwdAztecOOParamList.set("AZ_subdomain_solve",aztecSubdomainSolve);
        adjAztecOOParamList.set("AZ_subdomain_solve",aztecSubdomainSolve);
      }
    }
*/
    
    success
      = Thyra::test_single_aztecoo_thyra_solver(
        matrixFile,testTranspose,numRandomVectors
        ,maxFwdError,maxResid,maxSolutionError,showAllTests,dumpAll
        ,&paramList
        ,verbose?&out:0
        );

  }
  catch( const std::exception &excpt ) {
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    std::cerr << "*** Caught an unknown exception\n";
    success = false;
  }
  
  if (verbose) {
    if(success)  out << "\nCongratulations! All of the tests checked out!\n";
    else         out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? 0 : 1 );
}
