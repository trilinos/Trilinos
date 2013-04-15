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
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardCatchMacros.hpp>

int
main (int argc, char *argv[])
{
  using Teuchos::inOutArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::endl;
  typedef double                          ST;
  typedef Epetra_Operator                 OP;
  typedef Epetra_MultiVector              MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;

  // This calls MPI_Init and MPI_Finalize as necessary.
  Belos::Test::MPISession session (inOutArg (argc), inOutArg (argv));
  RCP<const Epetra_Comm> comm = session.getComm ();
  const int MyPID = comm->MyPID ();

  //
  // Parameters to read from command-line processor
  //
  bool verbose = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numRHS = 1;  // total number of right-hand sides to solve for
  int maxIters = 13000;  // maximum number of iterations for solver to use
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
  cmdp.setOption ("num-rhs", &numRHS, "Number of right-hand sides to solve.");
  cmdp.setOption ("max-iters", &maxIters, "Maximum number of iterations per "
      "linear system (-1 means \"adapt to problem/block size\").");
  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  Teuchos::oblackholestream blackHole;
  std::ostream& verbOut = (verbose && MyPID == 0) ? std::cout : blackHole;

  bool testPassed = false;
  try
  {
    //
    // Generate the linear system(s) to solve.
    //
    verbOut << "Generating the linear system(s) to solve" << endl << endl;
    RCP<Epetra_CrsMatrix> A;
    RCP<Epetra_MultiVector> B, X;
    RCP<Epetra_Map> rowMap;
    try {
      // This might change the number of right-hand sides, if we read in
      // a right-hand side from the Harwell-Boeing file.
      Belos::createEpetraProblem (comm, filename, rowMap, A, B, X, numRHS);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION (true, std::runtime_error,
          "Failed to create Epetra problem for matrix "
          "filename \"" << filename << "\".  "
          "createEpetraProblem() reports the following "
          "error: " << e.what());
    }
    //
    // Compute the initial residual norm of the problem, so we can see
    // by how much it improved after the solve.
    //
    std::vector<double> initialResidualNorms (numRHS);
    std::vector<double> initialResidualInfNorms (numRHS);
    Epetra_MultiVector R (*rowMap, numRHS);
    OPT::Apply (*A, *X, R);
    MVT::MvAddMv (-1.0, R, 1.0, *B, R); // R := -(A*X) + B.
    MVT::MvNorm (R, initialResidualNorms);
    MVT::MvNorm (R, initialResidualInfNorms, Belos::InfNorm);
    if (verbose) {
      verbOut << "Initial residual 2-norms:            \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << initialResidualNorms[i];
        if (i < numRHS-1) {
    verbOut << ", ";
        }
      }
      verbOut << endl << "Initial residual Inf-norms:          \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << initialResidualInfNorms[i];
        if (i < numRHS-1) {
    verbOut << ", ";
        }
      }
      verbOut << endl;
    }

    std::vector<double> rhs2Norms (numRHS);
    std::vector<double> rhsInfNorms (numRHS);
    MVT::MvNorm (*B, rhs2Norms);
    MVT::MvNorm (*B, rhsInfNorms, Belos::InfNorm);
    if (verbose) {
      verbOut << "Right-hand side 2-norms:             \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << rhs2Norms[i];
        if (i < numRHS-1) {
    verbOut << ", ";
        }
      }
      verbOut << endl << "Right-hand side Inf-norms:           \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << rhsInfNorms[i];
        if (i < numRHS-1) {
    verbOut << ", ";
        }
      }
      verbOut << endl;
    }

    std::vector<double> initialGuess2Norms (numRHS);
    std::vector<double> initialGuessInfNorms (numRHS);
    MVT::MvNorm (*X, initialGuess2Norms);
    MVT::MvNorm (*X, initialGuessInfNorms, Belos::InfNorm);
    if (verbose) {
      verbOut << "Initial guess 2-norms:               \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << initialGuess2Norms[i];
        if (i < numRHS-1) {
    verbOut << ", ";
        }
      }
      verbOut << endl << "Initial guess Inf-norms:             \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << initialGuessInfNorms[i];
        if (i < numRHS-1) {
    verbOut << ", ";
        }
      }
      verbOut << endl;
    }
    //
    // Compute the infinity-norm of A.
    //
    const double normOfA = A->NormInf ();
    verbOut << "||A||_inf:                           \t" << normOfA << endl;
    //
    // Compute ||A|| ||X_i|| + ||B_i|| for each right-hand side B_i.
    //
    std::vector<double> scaleFactors (numRHS);
    for (int i = 0; i < numRHS; ++i) {
      scaleFactors[i] = normOfA * initialGuessInfNorms[i] + rhsInfNorms[i];
    }
    if (verbose) {
      verbOut << "||A||_inf ||X_i||_inf + ||B_i||_inf: \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << scaleFactors[i];
        if (i < numRHS-1) {
    verbOut << ", ";
        }
      }
      verbOut << endl;
    }

    //
    // Solve using Belos
    //
    verbOut << endl << "Setting up Belos" << endl;
    const int NumGlobalElements = B->GlobalLength();

    // Set up Belos solver parameters.
    RCP<ParameterList> belosList = parameterList ("MINRES");
    belosList->set ("Maximum Iterations", maxIters);
    belosList->set ("Convergence Tolerance", tol);
    if (verbose) {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
          Belos::IterationDetails + Belos::OrthoDetails +
          Belos::FinalSummary + Belos::TimingDetails + Belos::Debug);
      belosList->set ("Output Frequency", frequency);
    }
    else {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings);
    }
    belosList->set ("Output Stream", rcpFromRef (verbOut));

    // Construct an unpreconditioned linear problem instance.
    typedef Belos::LinearProblem<double,MV,OP> prob_type;
    RCP<prob_type> problem = rcp (new prob_type (A, X, B));
    if (! problem->setProblem()) {
      verbOut << endl << "ERROR:  Failed to set up Belos::LinearProblem!" << endl;
      return EXIT_FAILURE;
    }

    // Create an iterative solver manager.
    RCP<Belos::SolverManager<double,MV,OP> > newSolver
      = rcp (new Belos::MinresSolMgr<double,MV,OP> (problem, belosList));

    // Print out information about problem.  Make sure to use the
    // information as stored in the Belos ParameterList, so that we know
    // what the solver will do.
    verbOut << endl
      << "Dimension of matrix: " << NumGlobalElements << endl
      << "Number of right-hand sides: " << numRHS << endl
      << "Max number of MINRES iterations: "
      << belosList->get<int> ("Maximum Iterations") << endl
      << "Relative residual tolerance: "
      << belosList->get<double> ("Convergence Tolerance") << endl
      << "Output frequency: "
      << belosList->get<int> ("Output Frequency") << endl
      << endl;

    // Solve the linear system.
    verbOut << "Solving the linear system" << endl << endl;
    Belos::ReturnType ret = newSolver->solve();
    verbOut << "Belos results:" << endl
      << "- Number of iterations: "
      << newSolver->getNumIters () << endl
      << "- " << (ret == Belos::Converged ? "Converged" : "Not converged")
      << endl;
    //
    // After the solve, compute residual(s) explicitly.  This tests
    // whether the Belos solver did so correctly.
    //
    std::vector<double> absoluteResidualNorms (numRHS);
    OPT::Apply (*A, *X, R);
    MVT::MvAddMv (-1.0, R, 1.0, *B, R);
    MVT::MvNorm (R, absoluteResidualNorms);

    std::vector<double> relativeResidualNorms (numRHS);
    for (int i = 0; i < numRHS; ++i) {
      relativeResidualNorms[i] = (initialResidualNorms[i] == 0.0) ?
        absoluteResidualNorms[i] :
        absoluteResidualNorms[i] / initialResidualNorms[i];
    }

    verbOut << "---------- Computed relative residual norms ----------"
      << endl << endl;
    bool badRes = false;
    if (verbose) {
      for (int i = 0; i < numRHS; ++i) {
        const double actRes = relativeResidualNorms[i];
        verbOut << "Problem " << i << " : \t" << actRes << endl;
        if (actRes > tol) {
    badRes = true;
        }
      }
    }

#   ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor::summarize (verbOut);
#   endif // BELOS_TEUCHOS_TIME_MONITOR

    testPassed = (ret == Belos::Converged && !badRes);
  } // try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, verbOut, testPassed);

  if (testPassed) {
    verbOut << endl << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  } else {
    verbOut << endl << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
} // end test_minres_hb.cpp
