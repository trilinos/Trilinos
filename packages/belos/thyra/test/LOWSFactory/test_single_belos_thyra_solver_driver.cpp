
#include "test_single_belos_thyra_solver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[])
{
  
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

	try {

    //
    // Read options from command-line
    //
    
    std::string     matrixFile             = "";
    bool            testTranspose          = true;
    int             numRandomVectors       = 1;
    double          maxFwdError            = 1e-14;
    int             maxIterations          = 400;
    int             maxRestarts            = 25;
    int             gmresKrylovLength      = 25;
    int             blockSize              = 1;
    double          maxResid               = 1e-6;
    double          maxSolutionError       = 1e-6;
    bool            showAllTests           = false;
    bool            dumpAll                = false;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "matrix-file", &matrixFile, "Matrix input file [Required]." );
    clp.setOption( "test-transpose", "no-test-transpose", &testTranspose, "Test the transpose solve or not." );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "max-fwd-error", &maxFwdError, "The maximum relative error in the forward operator." );
    clp.setOption( "max-iters", &maxIterations, "The maximum number of linear solver iterations to take." );
    clp.setOption( "max-restarts", &maxRestarts, "???." );
    clp.setOption( "gmres-krylov-length", &gmresKrylovLength, "???." );
    clp.setOption( "block-size", &blockSize, "???." );
    clp.setOption( "max-resid", &maxResid, "The maximum relative error in the residual." );
    clp.setOption( "max-solution-error", &maxSolutionError, "The maximum relative error in the solution of the linear system." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPT( matrixFile == "" );

    Teuchos::ParameterList solveParamList;

    solveParamList.set("Solver Type","GMRES");
    solveParamList.set("Max Iters",int(maxIterations));
    solveParamList.set("Max Restarts",int(maxRestarts));
    solveParamList.set("Block Size",int(blockSize));
    solveParamList.sublist("GMRES").set("Length",int(gmresKrylovLength));

    success
      = Thyra::test_single_belos_thyra_solver(
        matrixFile,testTranspose,numRandomVectors
        ,maxFwdError,maxResid,maxSolutionError,showAllTests,dumpAll
        ,&solveParamList
        ,verbose?&*out:0
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
		if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
		else         *out << "\nOh no! At least one of the tests failed!\n";
	}

  return ( success ? 0 : 1 );
}
