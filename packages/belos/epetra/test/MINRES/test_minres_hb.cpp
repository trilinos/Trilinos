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
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosMinresSolMgr.hpp"
#include "createEpetraProblem.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_MPI
#  include <mpi.h>
#endif // HAVE_MPI

int 
main (int argc, char *argv[]) 
{
  using Teuchos::inOutArg;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::endl;

  Belos::Test::MPISession session (inOutArg (argc), inOutArg (argv));
  RCP<const Epetra_Comm> comm = session.getComm ();
  const int MyPID = comm->MyPID ();

  //
  // Parameters to read from command-line processor
  // 
  bool verbose = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;  // total number of right-hand sides to solve for
  int blocksize = 1;  // blocksize used by solver
  int maxiters = -1;  // maximum number of iterations for solver to use
  std::string filename ("bcsstk14.hb");
  double tol = 1.0e-5; // relative residual tolerance
  
  //
  // Read in command-line arguments
  //
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("frequency", &frequency, "Solvers frequency for printing "
		  "residuals (#iters).");
  cmdp.setOption ("tol", &tol, "Relative residual tolerance used by MINRES "
		  "solver.");
  cmdp.setOption ("filename", &filename, "Filename for Harwell-Boeing test "
		  "matrix.");
  cmdp.setOption ("num-rhs", &numrhs, "Number of right-hand sides to solve.");
  cmdp.setOption ("max-iters", &maxiters, "Maximum number of iterations per "
		  "linear system (-1 means \"adapt to problem/block size\").");
  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (! verbose) {
    // In non-verbose mode, override the specified frequency so that
    // the solver doesn't print intermediate output.
    frequency = -1;  
  }
  //
  // Generate the linear system(s) to solve.
  //
  RCP<Epetra_CrsMatrix> A;
  RCP<Epetra_MultiVector> B, X;
  RCP<Epetra_Map> rowMap;
  try {
    Belos::createEpetraProblem (comm, filename, rowMap, A, B, X);
  } catch (std::exception& e) {
    TEUCHOS_TEST_FOR_EXCEPTION (true, std::runtime_error,
				"Failed to create Epetra problem for matrix "
				"filename \"" << filename << "\".  "
				"createEpetraProblem() reports the following "
				"error: " << e.what());
  }

  Teuchos::oblackholestream blackHole;
  std::ostream& verbOut = (verbose && MyPID == 0) ? std::cout : blackHole;

  //
  // Solve using Belos
  //
  typedef double                          ST;
  typedef Epetra_Operator                 OP;
  typedef Epetra_MultiVector              MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  //
  // *****Construct initial guess and random right-hand-sides *****
  //
  if (numrhs != 1) {
    X = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
    MVT::MvRandom( *X );
    B = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );
  }
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  //
  ParameterList belosList;
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  if (verbose) {
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + 
		   Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  else
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

  // Construct an unpreconditioned linear problem instance.
  Belos::LinearProblem<double,MV,OP> problem (A, X, B);
  if (! problem.setProblem()) {
    verbOut << endl << "ERROR:  Failed to set up Belos::LinearProblem!" << endl;
    return EXIT_FAILURE;
  }

  // Create an iterative solver manager.
  RCP<Belos::SolverManager<double,MV,OP> > newSolver
    = rcp (new Belos::MinresSolMgr<double,MV,OP> (rcpFromRef (problem), 
						  rcpFromRef (belosList)));
  // Print out information about problem.
  verbOut << endl << endl
	  << "Dimension of matrix: " << NumGlobalElements << endl
	  << "Number of right-hand sides: " << numrhs << endl
	  << "Max number of MINRES iterations: " << maxiters << endl
	  << "Relative residual tolerance: " << tol << endl
	  << endl;
  // Solve the linear system.
  Belos::ReturnType ret = newSolver->solve();
  //
  // Compute residual(s) explicitly.  This tests whether the Belos
  // solver did so correctly.
  //
  bool badRes = false;
  std::vector<double> actual_resids (numrhs);
  std::vector<double> rhs_norm (numrhs);
  Epetra_MultiVector resid (A->Map(), numrhs);
  OPT::Apply (*A, *X, resid);
  MVT::MvAddMv (-1.0, resid, 1.0, *B, resid); 
  MVT::MvNorm (resid, actual_resids);
  MVT::MvNorm (*B, rhs_norm);

  verbOut << "---------- Actual Residuals (normalized) ----------" 
	  << endl << endl;
  if (verbose) {
    for (int i = 0; i < numrhs; ++i) {
      double actRes = actual_resids[i]/rhs_norm[i];
      verbOut << "Problem " << i << " : \t" << actRes << endl;
      if (actRes > tol) {
	badRes = true;
      }
    }
  }

  const bool testFailed = (ret != Belos::Converged || badRes);
  if (testFailed) {
    verbOut << endl << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  } else {
    verbOut << endl << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
} // end test_minres_hb.cpp




