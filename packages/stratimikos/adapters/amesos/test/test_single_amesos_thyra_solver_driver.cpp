/*
// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "test_single_amesos_thyra_solver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[])
{
  
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //
    
    std::string                             matrixFile             = "";
    Thyra::Amesos::ESolverType              solverType
#ifdef HAVE_AMESOS_KLU
                                                                   = Thyra::Amesos::KLU;
#else
                                                                   = Thyra::Amesos::LAPACK;
#endif
    Thyra::Amesos::ERefactorizationPolicy   refactorizationPolicy  = Thyra::Amesos::REPIVOT_ON_REFACTORIZATION;
    bool                                    testTranspose          = true;
    int                                     numRandomVectors       = 1;
    double                                  maxFwdError            = 1e-14;
    double                                  maxError               = 1e-10;
    double                                  maxResid               = 1e-10;
    bool                                    showAllTests           = false;
    bool                                    dumpAll                = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "matrix-file", &matrixFile, "Matrix iput file [Required]." );
    clp.setOption(
      "solver-type", &solverType
      ,Thyra::Amesos::numSolverTypes, Thyra::Amesos::solverTypeValues, Thyra::Amesos::solverTypeNames
      ,"Type of direct solver."
      );
    clp.setOption(
      "refactorization-policy", &refactorizationPolicy
      ,Thyra::Amesos::numRefactorizationPolices, Thyra::Amesos::refactorizationPolicyValues, Thyra::Amesos::refactorizationPolicyNames
      ,"Pivoting policy used on refactorizations."
      );
    clp.setOption( "test-transpose", "no-test-transpose", &testTranspose, "Test the transpose solve or not." );
    clp.setOption( "num-random-vectors", &numRandomVectors, "Number of times a test is performed with different random vectors." );
    clp.setOption( "max-fwd-error", &maxFwdError, "The maximum relative error in the forward operator." );
    clp.setOption( "max-error", &maxError, "The maximum relative error in the solution." );
    clp.setOption( "max-resid", &maxResid, "The maximum relative error in the residual." );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Set if all the tests are shown or not." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEUCHOS_TEST_FOR_EXCEPT( matrixFile == "" );

    Teuchos::ParameterList amesosLOWSFPL;
    amesosLOWSFPL.set("Solver Type",toString(solverType));
    amesosLOWSFPL.set("Refactorization Policy",toString(refactorizationPolicy));

    success
      = Thyra::test_single_amesos_thyra_solver(
        matrixFile,&amesosLOWSFPL,testTranspose,numRandomVectors
        ,maxFwdError,maxError,maxResid,showAllTests,dumpAll,verbose?&*out:0
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
