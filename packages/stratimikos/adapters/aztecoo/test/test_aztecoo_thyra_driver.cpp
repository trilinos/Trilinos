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
#include "Teuchos_StandardCatchMacros.hpp"
#include "az_aztec_defs.h"

struct MatrixTestPacket {
  MatrixTestPacket(
    std::string  _matrixFile
    ,double      _maxFwdError
    ,int         _maxIters
    ,double      _maxResid
    ,double      _maxSolutionError
    ,double      _maxSlackErrorFrac
    ,int         _maxPrecIters
    ,double      _maxPrecResid
    ,double      _maxPrecSolutionError
    ,double      _maxPrecSlackErrorFrac
    )
    :matrixFile(_matrixFile)
    ,maxFwdError(_maxFwdError)
    ,maxIters(_maxIters)
    ,maxResid(_maxResid)
    ,maxSolutionError(_maxSolutionError)
    ,maxSlackErrorFrac(_maxSlackErrorFrac)
    ,maxPrecIters(_maxPrecIters)
    ,maxPrecResid(_maxPrecResid)
    ,maxPrecSolutionError(_maxPrecSolutionError)
    ,maxPrecSlackErrorFrac(_maxPrecSlackErrorFrac)
    {}
  std::string  matrixFile;
  double       maxFwdError;
  int          maxIters;
  double       maxResid;
  double       maxSolutionError;
  double       maxSlackErrorFrac;
  int          maxPrecIters;
  double       maxPrecResid;
  double       maxPrecSolutionError;
  double       maxPrecSlackErrorFrac;
};

int main(int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  using Teuchos::CommandLineProcessor;

  bool result, success = true;
  bool verbose = true;

  Teuchos::FancyOStream out(Teuchos::rcp(&std::cout,false));

  try {

    //
    // Read options from command-line
    //

    std::string    matrixDir              = ".";
    int            numRandomVectors       = 1;
    bool           showAllTests           = false;
    bool           showAllTestsDetails    = false;
    bool           dumpAll                = false;
    std::string    aztecOutputLevel       = "freq";
    int            aztecOutputFreq        = 0;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "matrix-dir", &matrixDir, "Base directory for the test matrices" );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "show-all-tests-details", "no-show-all-tests-details", &showAllTestsDetails, "Set if all the details of the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    clp.setOption( "aztec-output-level", &aztecOutputLevel, "Aztec output level (freq,last,summary,warnings,all)" );
    clp.setOption( "aztec-output-freq", &aztecOutputFreq, "Aztec output freqency (> 0)" );
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
        MTP("bcsstk01.mtx"       ,1e-12, 40 , 1e-4, 0.6,      1.0, 20 , 1e-10, 0.5,      1.0)
        ,MTP("bcsstk02.mtx"      ,1e-12, 40 , 1e-3, 0.5,      1.0, 2  , 1e-10, 0.5,      1.0)
        ,MTP("bcsstk04.mtx"      ,1e-12, 80 , 1e-4, 0.999990, 1.0, 40 , 1e-10, 0.999990, 1.0)
        ,MTP("Diagonal.mtx"      ,1e-12, 4  , 1e-6, 1e-14,    1.0, 2  , 1e-10, 1e-14,    1.0)
        ,MTP("FourByFour.mtx"    ,1e-12, 4  , 1e-6, 1e-14,    1.0, 2  , 1e-10, 1e-14,    1.0)
        ,MTP("KheadK.mtx"        ,1e-12, 8  , 1e-6, 1e-14,    1.0, 2  , 1e-10, 1e-14,    1.0)
        ,MTP("KheadSorted.mtx"   ,1e-12, 8  , 1e-6, 1e-14,    1.0, 2  , 1e-10, 1e-14,    1.0)
        ,MTP("nos1.mtx"          ,1e-11, 200, 1e-4, 0.8,      1.0, 237, 1e-2,  5.0,      1.0)
        ,MTP("nos5.mtx"          ,1e-12, 468, 1e-5, 0.5,      1.0, 468, 1e-10, 0.5,      1.0)
      };
    //
    // Loop through all of the test matrices
    //
    for( int matrix_i = 0; matrix_i < numTestMatrices; ++matrix_i ) {
      const MatrixTestPacket
        mtp = testMatrices[matrix_i];
      //
      // Do unpreconditioned and preconditioned solves
      //
      for( int prec_i = 0; prec_i < 2; ++prec_i ) {
        if(verbose)
          out << std::endl<<matrix_i<<":"<<prec_i<<": Testing, matrixFile=\'"<<mtp.matrixFile<<"\', ";
        bool testTranspose;
        double maxResid;
        double maxSolutionError;
        //double maxSlackErrorFrac;
        Teuchos::ParameterList
          paramList("AztecOOLinearOpWithSolveFactory");
        Teuchos::ParameterList
          &fwdSolvePL = paramList.sublist("Forward Solve"),
          &adjSolvePL = paramList.sublist("Adjoint Solve");
        Teuchos::ParameterList
          &fwdAztecOOPL = fwdSolvePL.sublist("AztecOO Settings"),
          &adjAztecOOPL = adjSolvePL.sublist("AztecOO Settings");
        if( aztecOutputLevel != "freq" ) {
          fwdAztecOOPL.set("Output Frequency",aztecOutputLevel);
          adjAztecOOPL.set("Output Frequency",aztecOutputLevel);
        }
        else {
          fwdAztecOOPL.set("Output Frequency",aztecOutputFreq);
          adjAztecOOPL.set("Output Frequency",aztecOutputFreq);
        }
        if(prec_i==0) {
          out << "no aztec preconditioning ... ";
          fwdAztecOOPL.set("Aztec Preconditioner","none");
          testTranspose = true;
          fwdSolvePL.set("Max Iterations",mtp.maxIters);
          adjSolvePL.set("Max Iterations",mtp.maxIters);
          maxResid = mtp.maxResid;
          maxSolutionError = mtp.maxSolutionError;
          //maxSlackErrorFrac = mtp.maxSlackErrorFrac;
        }
        else {
          out << "using aztec preconditioning ... ";
          fwdAztecOOPL.set("Aztec Preconditioner","ilu");
          testTranspose = false;
          fwdSolvePL.set("Max Iterations",mtp.maxPrecIters);
          adjSolvePL.set("Max Iterations",mtp.maxPrecIters);
          maxResid = mtp.maxPrecResid;
          maxSolutionError = mtp.maxPrecSolutionError;
          //maxSlackErrorFrac = mtp.maxPrecSlackErrorFrac;
        }
        std::ostringstream oss;
        Teuchos::FancyOStream fancy_oss(Teuchos::rcp(&oss,false));
        result =
          Thyra::test_single_aztecoo_thyra_solver(
            matrixDir+"/"+mtp.matrixFile,testTranspose,numRandomVectors
            ,mtp.maxFwdError,maxResid,maxSolutionError
            ,showAllTestsDetails,dumpAll,&paramList,&fancy_oss
            );
        if(!result) success = false;
        if(verbose) {
          if(result) {
            if(showAllTests)
              out << std::endl << oss.str();
            else
              out << " : passed!\n";
          }
          else {
            if(showAllTests)
              out << std::endl << oss.str();
            else
              out << " : failed!\n";
          }
        }
      }
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (verbose) {
    if(success)  out << "\nCongratulations! All of the tests checked out!\n";
    else         out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
