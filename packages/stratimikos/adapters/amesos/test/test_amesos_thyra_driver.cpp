// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_single_amesos_thyra_solver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"

struct MatrixTestPacket {
  MatrixTestPacket(std::string  _matrixFile, bool _unsymmetric, double _maxFwdError, double _maxError, double _maxResid)
    :matrixFile(_matrixFile),unsymmetric(_unsymmetric),maxFwdError(_maxFwdError),maxError(_maxError),maxResid(_maxResid) {}
  std::string  matrixFile;
  bool         unsymmetric;
  double       maxFwdError;
  double       maxError;
  double       maxResid;
};

int main(int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  using Teuchos::CommandLineProcessor;
  using Teuchos::OSTab;

  bool result, success = true;
  bool verbose = true;


  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //
    
    std::string    matrixDir              = ".";
    bool           testTranspose          = true;
    int            numRandomVectors       = 1;
    bool           showAllTests           = false;
    bool           showAllTestsDetails    = false;
    bool           dumpAll                = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "matrix-dir", &matrixDir, "Base directory for the test matrices" );
    clp.setOption( "test-transpose", "no-test-transpose", &testTranspose, "Test the transpose solve or not." );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "show-all-tests-details", "no-show-all-tests-details", &showAllTestsDetails, "Set if all the details of the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEUCHOS_TEST_FOR_EXCEPT( matrixDir == "" );

    //
    // Define the test matrices
    //

    const int numTestMatrices = 9;

    typedef MatrixTestPacket MTP;

    // Set up the matices and the tolerances.
    // Note, we may need to adjust these for bad platforms ...
    const MTP testMatrices[numTestMatrices] =
      {
        MTP("bcsstk01.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("bcsstk02.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("bcsstk04.mtx",false,1e-12,1e-10,1e-12)
        ,MTP("Diagonal.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("FourByFour.mtx",true,1e-12,1e-12,1e-12)
        ,MTP("KheadK.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("KheadSorted.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("nos1.mtx",false,1e-11,1e-10,1e-12)
        ,MTP("nos5.mtx",false,1e-12,1e-12,1e-12)
      };
    //
    // Loop through all of the test matrices
    //
    for( int matrix_i = 0; matrix_i < numTestMatrices; ++matrix_i ) {
      const MatrixTestPacket
        mtp = testMatrices[matrix_i];
      //
      // Loop through all of the solvers
      //
      for( int solver_i = 0; solver_i < Thyra::Amesos::numSolverTypes; ++solver_i ) {
        const Thyra::Amesos::ESolverType
          solverType = Thyra::Amesos::solverTypeValues[solver_i];

        //  bug 1902 - Amesos_Superlu fails on bcsstk01.mtx
        //  bug 1903 - Amesos_Superlu fails on four matrices,
        //             when called from the thyra test
        //
        bool BadMatrixForSuperlu = 
          mtp.matrixFile == "bcsstk01.mtx" // bug 1902 
          || mtp.matrixFile == "bcsstk04.mtx" // bug 1903 
          || mtp.matrixFile == "KheadK.mtx" // bug 1903 
          || mtp.matrixFile == "KheadSorted.mtx" // bug 1903 
          || mtp.matrixFile == "nos1.mtx" ; // bug 1903 
        //
        // Toggle the refactorization options
        //
        for( int factorizationPolicy_i = 0; factorizationPolicy_i < Thyra::Amesos::numRefactorizationPolices;  ++factorizationPolicy_i ) {
          const Thyra::Amesos::ERefactorizationPolicy
            refactorizationPolicy = Thyra::Amesos::refactorizationPolicyValues[factorizationPolicy_i];
          if(verbose)
            *out
              << std::endl<<matrix_i<<"."<<solver_i<<"."<<factorizationPolicy_i<<": "
              << "Testing, matrixFile=\'"<<mtp.matrixFile<<"\', solverType=\'"<<toString(solverType)<<"\', refactorizationPolicy=\'"<<toString(refactorizationPolicy)<<"\' ..."; 
          if( mtp.unsymmetric && !Thyra::Amesos::supportsUnsymmetric[solver_i] ) {
            *out << " : Skipping since unsymmetric and not supported!\n";
          }
          else {
            //  bug 1902 and bug 1903  
            std::string StrSolverType = toString(solverType) ; 
            std::string StrSuperlu = "Superlu";
            if ( StrSolverType==StrSuperlu && BadMatrixForSuperlu ) {
              *out << " : Skipping since Superlu fails on this matrix!\n";
            }
            else {
              std::ostringstream ossStore;
              Teuchos::RCP<Teuchos::FancyOStream>
                oss = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&ossStore,false)));
              Teuchos::ParameterList amesosLOWSFPL;
              amesosLOWSFPL.set("Solver Type",toString(solverType));
              amesosLOWSFPL.set("Refactorization Policy",toString(refactorizationPolicy));
              result =
                Thyra::test_single_amesos_thyra_solver(
                  matrixDir+"/"+mtp.matrixFile,&amesosLOWSFPL,testTranspose,numRandomVectors
                  ,mtp.maxFwdError,mtp.maxError,mtp.maxResid,showAllTestsDetails,dumpAll,OSTab(oss).get()
                  );
              if(!result) success = false;
              if(verbose) {
                if(result) {
                  if(showAllTests)
                    *out << std::endl << ossStore.str();
                  else
                    *out << " : passed!\n";
                }
                else {
                  if(showAllTests)
                    *out << std::endl << ossStore.str();
                  else
                    *out << " : failed!\n";
                }
              }
            }
          }
        }
      }
    }
    
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
