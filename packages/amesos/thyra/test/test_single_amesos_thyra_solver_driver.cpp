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

    TEST_FOR_EXCEPT( matrixFile == "" );

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
