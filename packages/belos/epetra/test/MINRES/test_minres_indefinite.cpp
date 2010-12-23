//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
//
// This driver reads a problem from a file, which can be in Harwell-Boeing (*.hb),
// Matrix Market (*.mtx), or triplet format (*.triU, *.triS).  The right-hand side
// from the problem, if it exists, will be used instead of multiple
// random right-hand-sides.  The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this example. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosMinresSolMgr.hpp"

#include "Epetra_Map.h"
#ifdef EPETRA_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

int 
main(int argc, char *argv[]) 
{
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;

  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  const int MyPID = Comm.MyPID();
#else // No EPETRA_MPI
  Epetra_SerialComm Comm;
  const int MyPID = 0;
#endif // EPETRA_MPI

  bool verbose = false, debug = false;
  int frequency = -1;        // frequency of status test output.
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  MT tol = 1.0e-8;           // relative residual tolerance

  // Define command-line arguments
  CommandLineProcessor cmdp(false,true);
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug, 
		  "Print debugging information from the solver.");
  cmdp.setOption ("frequency", &frequency,
		  "Solver's frequency for printing residuals (#iters).");
  cmdp.setOption ("tol", &tol,
		  "Relative residual tolerance used by MINRES solver.");
  cmdp.setOption ("num-rhs", &numrhs,
		  "Number of right-hand sides for which to solve (> 0).");
  cmdp.setOption ("max-iters", &maxiters, 
		  "Maximum number of iterations per linear system.  -1 means "
		  "we choose, based on problem size.");

  // Parse command-line arguments and fetch values
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) 
    {
      std::cout << "Failed to parse command-line arguments!" << std::endl;
      std::cout << "End Result: TEST FAILED" << std::endl;
      return -1;
    }

  // **********************************************************************
  // ******************Set up the problem to be solved*********************
  // construct diagonal matrix
  const int NumGlobalElements = 1e2;
  const int m = 4; // number of negative eigenvalues
  
  // Create diagonal matrix with n-m positive and m negative eigenvalues.
  Epetra_Map epetraMap( NumGlobalElements, 0, Comm );
  
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix( Copy, epetraMap, 1 ) );
  for ( int k=0; k<epetraMap.NumMyElements(); k++ )
  {
      int GIDk = epetraMap.GID(k);
      double val = 2*(GIDk-m) + 1;
      TEUCHOS_ASSERT_EQUALITY( 0, A->InsertGlobalValues( GIDk, 1, &val, &GIDk ) );
  } 
  TEUCHOS_ASSERT_EQUALITY( 0, A->FillComplete() );
  TEUCHOS_ASSERT_EQUALITY( 0, A->OptimizeStorage() );
  
  //
  // Make some (multi)vectors, for use in testing: an exact solution,
  // its corresponding right-hand side, and an initial guess. 
  //

  // Make a random exact solution.
  Teuchos::RCP<Epetra_MultiVector> X_exact (new Epetra_MultiVector (epetraMap, numrhs));
  X_exact->Seed();
  MVT::MvRandom (*X_exact);

  // Compute the right-hand side as B = A*X.
  Teuchos::RCP<Epetra_MultiVector> B = MVT::Clone (*X_exact, numrhs);
  OPT::Apply (*A, *X_exact, *B);

  // Choose an initial guess of all zeros.
  Teuchos::RCP<Epetra_MultiVector> X = MVT::Clone (*X_exact, numrhs);
  MVT::MvInit (*X, ST(0.0));

  // **********************************************************************
  //
  // Set parameters that the user asked us to pick intelligently.
  //

  // Maximum number of iterations: set according to problem size.
  // maxiters=numGlobalElements-1 may not always converge, so we
  // set to numGlobalElements+1 for good measure.
  if (maxiters == -1)
    maxiters = NumGlobalElements + 1; 

  // In a nonverbose test, the frequency should be set to -1, which
  // Belos interprets as "no intermediate status output."  Override
  // whatever the user may have specified.
  if (! verbose)
    frequency = -1;
  // Silently fix a bad frequency value.
  else if (frequency < 0 && frequency != -1)
    frequency = -1;

  // Validate command-line arguments
  TEST_FOR_EXCEPTION( tol < 0, std::invalid_argument,
		      "Relative residual tolerance must be nonnegative, but "
		      "you supplied tol = " << tol << "." );
  TEST_FOR_EXCEPTION( numrhs < 1, std::invalid_argument,
		      "MINRES test requires at least one right-hand side, but "
		      "you set the number of right-hand sides to " 
		      << numrhs << "." );
  TEST_FOR_EXCEPTION( maxiters < 1, std::invalid_argument,
		      "MINRES test requires at least one iteration, but you "
		      "set the maximum number of iterations to " 
		      << maxiters << "." );

  const bool proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //

  //
  ParameterList belosList;
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested

  int verbosity = Belos::Errors + Belos::Warnings;
  if (verbose) {
    verbosity += Belos::TimingDetails + Belos::StatusTestDetails;
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  if (debug) {
    verbosity += Belos::Debug;
  }
  belosList.set( "Verbosity", (int) verbosity );
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem< ST, MV, OP > problem (A, X, B);
  {
    const bool set = problem.setProblem();
    TEST_FOR_EXCEPTION( set == false, std::logic_error, 
			"Belos::LinearProblem failed to set up correctly (setP"
			"roblem() returned false)!  This probably means we imp"
			"lemented our test incorrectly." );
  }
  //
  // *******************************************************************
  // ****************Start the MINRES iteration*************************
  // *******************************************************************
  //
  Belos::OutputManager< ST > My_OM();
 
  // Create an iterative solver manager.
  RCP< Belos::SolverManager<ST,MV,OP> > newSolver
    = rcp( new Belos::MinresSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    cout << endl << endl;
    cout << "Dimension of matrix: " << NumGlobalElements << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }
  //
  // Perform solve
  //
  Belos::ReturnType ret = Belos::Unconverged;
  try {
    ret = newSolver->solve();
  } catch (std::exception& e) {
    if (MyPID == 0)
      cout << "MINRES solver threw an exception: " << e.what() << endl;
    throw e;
  }
  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(A->Map(), numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    for ( int i=0; i<numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      if (actRes > tol) badRes = true;
    }
  }

  if (ret!=Belos::Converged || badRes) {
    if (proc_verbose)
      std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
  }
  //
  // Default return value
  //
  if (proc_verbose)
    std::cout << std::endl << "End Result: TEST PASSED" << std::endl;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  if (ret!=Belos::Converged || badRes)
    return -1;
  else
    return 0;
} 
