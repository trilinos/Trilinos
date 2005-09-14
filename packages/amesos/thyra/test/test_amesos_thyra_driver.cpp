/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "test_single_amesos_thyra_solver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

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
  
  using Teuchos::CommandLineProcessor;

  bool result, success = true;
  bool verbose = true;

  std::ostream &out = std::cout;

	try {

    //
    // Read options from command-line
    //
    
    std::string    matrixDir              = "";
    bool           testTranspose          = true;
    int            numRandomVectors       = 1;
    bool           showAllTests           = false;
    bool           showAllTestsDetails    = false;
    bool           dumpAll                = false;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "matrix-dir", &matrixDir, "Base directory for the test matrices" );
    clp.setOption( "test-transpose", "no-test-transpose", &testTranspose, "Test the transpose solve or not." );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "show-all-tests-details", "no-show-all-tests-details", &showAllTestsDetails, "Set if all the details of the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPT( matrixDir == "" );

    //
    // Define the test matrices
    //

    const int numTestMatrices = 9;

    typedef MatrixTestPacket MTP;

    // Set up the matices and the tolerances.
    // Note, we may need to adjust these for bad platforms ...
    const MTP testMatrices[numTestMatrices] =
      {
        MTP("In_bcsstk01.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("In_bcsstk02.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("In_bcsstk04.mtx",false,1e-12,1e-10,1e-12)
        ,MTP("In_Diagonal.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("In_FourByFour.mtx",true,1e-12,1e-12,1e-12)
        ,MTP("In_KheadK.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("In_KheadSorted.mtx",false,1e-12,1e-12,1e-12)
        ,MTP("In_nos1.mtx",false,1e-11,1e-10,1e-12)
        ,MTP("In_nos5.mtx",false,1e-12,1e-12,1e-12)
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
        //
        // Toggle the refactorization options
        //
        for( int factorizationPolicy_i = 0; factorizationPolicy_i < Thyra::Amesos::numRefactorizationPolices;  ++factorizationPolicy_i ) {
          const Thyra::Amesos::ERefactorizationPolicy
            refactorizationPolicy = Thyra::Amesos::refactorizationPolicyValues[factorizationPolicy_i];
          if(verbose)
            out << std::endl<<matrix_i<<"."<<solver_i<<"."<<factorizationPolicy_i<<": "
                << "Testing, matrixFile=\'"<<mtp.matrixFile<<"\', solverType=\'"<<toString(solverType)<<"\', refactorizationPolicy=\'"<<toString(refactorizationPolicy)<<"\' ..."; 
          if( mtp.unsymmetric && !Thyra::Amesos::supportsUnsymmetric[solver_i] ) {
            out << " : Skipping since unsymmetric and not supported!\n";
          }
          else {
            std::ostringstream oss;
            result =
              Thyra::test_single_amesos_thyra_solver(
                matrixDir+"/"+mtp.matrixFile,solverType,refactorizationPolicy,testTranspose,numRandomVectors
                ,mtp.maxFwdError,mtp.maxError,mtp.maxResid,showAllTestsDetails,dumpAll,&oss
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
		if(success)  out << "\nCongratulations! All of the tests checked out!\n";
		else         out << "\nOh no! At least one of the tests failed!\n";
	}

  return ( success ? 0 : 1 );
}
