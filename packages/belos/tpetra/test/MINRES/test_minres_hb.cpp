// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* Originally convert test here: belos/epetra/test/MINRES/test_minres_hb.cpp */

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MatrixIO.hpp> // I/O for Harwell-Boeing files

// Teuchos
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"

template<class ScalarType>
int run(int argc, char *argv[]) {
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;

  using OP = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;

  using tcrsmatrix_t = typename Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using Teuchos::inOutArg;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  // This calls MPI_Init and MPI_Finalize as necessary.
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  bool success = false;
  bool verbose = false;

  try {
    int MyPID = rank(*comm);

    //
    // Parameters to read from command-line processor
    //
    int frequency = -1;  // how often residuals are printed by solver
    int numRHS = 1;  // total number of right-hand sides to solve for
    int maxIters = 13000;  // maximum number of iterations for solver to use
    std::string filename ("bcsstk14.hb");
    ST tol = sqrt(std::numeric_limits<ST>::epsilon()); // relative residual tolerance
    
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
      return EXIT_FAILURE;
    }
    Teuchos::oblackholestream blackHole;
    std::ostream& verbOut = (verbose && MyPID == 0) ? std::cout : blackHole;
    
    //
    // Generate the linear system(s) to solve.
    //
    verbOut << "Generating the linear system(s) to solve" << std::endl << std::endl;    
    RCP<tcrsmatrix_t> A;
    Tpetra::Utils::readHBMatrix(filename, comm, A); 
    RCP<const Tpetra::Map<> > rowMap = A->getDomainMap();
  
    //
    // *****Construct initial guess and random right-hand-sides *****
    //
    RCP<MV> B, X;    
    X = rcp( new MV(rowMap, numRHS) );
    MVT::MvRandom( *X );
    B = rcp( new MV(rowMap, numRHS ) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );
    
    //
    // Compute the initial residual norm of the problem, so we can see
    // by how much it improved after the solve.
    //
    std::vector<ST> initialResidualNorms (numRHS);
    std::vector<ST> initialResidualInfNorms (numRHS);
    MV R(rowMap, numRHS);
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
      verbOut << std::endl << "Initial residual Inf-norms:          \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << initialResidualInfNorms[i];
        if (i < numRHS-1) {
          verbOut << ", ";
        }
      }
      verbOut << std::endl;
    }
    
    std::vector<ST> rhs2Norms (numRHS);
    std::vector<ST> rhsInfNorms (numRHS);
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
      verbOut << std::endl << "Right-hand side Inf-norms:           \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << rhsInfNorms[i];
        if (i < numRHS-1) {
          verbOut << ", ";
        }
      }
      verbOut << std::endl;
    }
    
    std::vector<ST> initialGuess2Norms (numRHS);
    std::vector<ST> initialGuessInfNorms (numRHS);
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
      verbOut << std::endl << "Initial guess Inf-norms:             \t";
      for (int i = 0; i < numRHS; ++i) {
        verbOut << initialGuessInfNorms[i];
        if (i < numRHS-1) {
          verbOut << ", ";
        }
      }
      verbOut << std::endl;
    }
    
    //
    // Compute the infinity-norm of A.
    //
    const ST normOfA = A->getFrobeniusNorm(); 
    verbOut << "||A||_inf:                           \t" << normOfA << std::endl;
    //
    // Compute ||A|| ||X_i|| + ||B_i|| for each right-hand side B_i.
    //
    std::vector<ST> scaleFactors (numRHS);
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
      verbOut << std::endl;
    }
    
    //
    // Solve using Belos
    //
    verbOut << std::endl << "Setting up Belos" << std::endl;
    const int NumGlobalElements = B->getGlobalLength();

    // Set up Belos solver parameters.
    RCP<ParameterList> belosList = Teuchos::parameterList ("MINRES");
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
    belosList->set ("Output Stream", Teuchos::rcpFromRef (verbOut));

    // Construct an unpreconditioned linear problem instance.
    typedef Belos::LinearProblem<ST,MV,OP> prob_type;
    RCP<prob_type> problem = rcp (new prob_type (A, X, B));
    if (! problem->setProblem()) {
      verbOut << std::endl << "ERROR:  Failed to set up Belos::LinearProblem!" << std::endl;
      return EXIT_FAILURE;
    }
    
    // Create an iterative solver manager.
    Belos::SolverFactory<ST, MV, OP> factory;
    RCP<Belos::SolverManager<ST,MV,OP> > newSolver =
      factory.create ("MINRES", belosList);
    newSolver->setProblem (problem);

    // Print out information about problem.  Make sure to use the
    // information as stored in the Belos ParameterList, so that we know
    // what the solver will do.
    verbOut << std::endl
      << "Dimension of matrix: " << NumGlobalElements << std::endl
      << "Number of right-hand sides: " << numRHS << std::endl
      << "Max number of MINRES iterations: "
      << belosList->get<int> ("Maximum Iterations") << std::endl
      << "Relative residual tolerance: "
      << belosList->get<double> ("Convergence Tolerance") << std::endl
      << "Output frequency: "
      << belosList->get<int> ("Output Frequency") << std::endl
      << std::endl;

    // Solve the linear system.
    verbOut << "Solving the linear system" << std::endl << std::endl;
    Belos::ReturnType ret = newSolver->solve();
    verbOut << "Belos results:" << std::endl
      << "- Number of iterations: "
      << newSolver->getNumIters () << std::endl
      << "- " << (ret == Belos::Converged ? "Converged" : "Not converged")
      << std::endl;
    
    //
    // After the solve, compute residual(s) explicitly.  This tests
    // whether the Belos solver did so correctly.
    //
    std::vector<ST> absoluteResidualNorms (numRHS);
    OPT::Apply (*A, *X, R);
    MVT::MvAddMv (-1.0, R, 1.0, *B, R);
    MVT::MvNorm (R, absoluteResidualNorms);

    std::vector<ST> relativeResidualNorms (numRHS);
    for (int i = 0; i < numRHS; ++i) {
      relativeResidualNorms[i] = (initialResidualNorms[i] == 0.0) ?
        absoluteResidualNorms[i] :
        absoluteResidualNorms[i] / initialResidualNorms[i];
    }

    verbOut << "---------- Computed relative residual norms ----------"
      << std::endl << std::endl;
    bool badRes = false;
    if (verbose) {
      for (int i = 0; i < numRHS; ++i) {
        const double actRes = relativeResidualNorms[i];
        verbOut << "Problem " << i << " : \t" << actRes << std::endl;
        if (actRes > tol) {
          badRes = true;
        }
      }
    }
    success = (ret == Belos::Converged && !badRes);
    
    if (success) {
      verbOut << std::endl << "End Result: TEST PASSED" << std::endl;
    } else {
      verbOut << std::endl << "End Result: TEST FAILED" << std::endl;
    }
  } // try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
} // end test_minres_hb.cpp

int main(int argc, char *argv[]) {
  // run with different scalar types
  return run<double>(argc, argv);
  // return run<float>(argc, argv); // FAILS
}
=======
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosMinresSolMgr.hpp"

// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

using namespace Teuchos;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using std::endl;
using std::cout;
using std::vector;
using Teuchos::tuple;

int main(int argc, char *argv[]) {
  typedef Tpetra::MultiVector<>::scalar_type ST;
  typedef ScalarTraits<ST>                SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::Operator<ST>             OP;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;

  GlobalMPISession mpisess(&argc,&argv,&cout);

  bool success = false;
  bool verbose = false;
  try {
    const ST one  = SCT::one();

    int MyPID = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool debug = false;
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;      // total number of right-hand sides to solve for
    int blocksize = 1;   // blocksize used by solver
    int maxiters = -1;   // maximum number of iterations for solver to use
    std::string filename("bcsstk14.hb");
    MT tol = 1.0e-5;     // relative residual tolerance

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
    cmdp.setOption("block-size",&blocksize,"Block size to be used by the CG solver.");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (debug) {
      verbose = true;
    }
    if (!verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }

    MyPID = rank(*comm);
    proc_verbose = ( verbose && (MyPID==0) );

    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }
    //
    // Get the data from the HB file and build the Map,Matrix
    //
    RCP<CrsMatrix<ST> > A;
    Tpetra::Utils::readHBMatrix(filename,comm,A);
    RCP<const Tpetra::Map<> > map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp( new MV(map,numrhs) );
    MVT::MvRandom( *X );
    B = rcp( new MV(map,numrhs) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );

    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1) {
      maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
    }
    //
    ParameterList belosList;
    belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    int verbLevel = Belos::Errors + Belos::Warnings;
    if (debug) {
      verbLevel += Belos::Debug;
    }
    if (verbose) {
      verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
    }
    belosList.set( "Verbosity", verbLevel );
    if (verbose) {
      if (frequency > 0) {
        belosList.set( "Output Frequency", frequency );
      }
    }
    //
    // Construct an unpreconditioned linear problem instance.
    //
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // *******************************************************************
    // *************Start the block CG iteration***********************
    // *******************************************************************
    //
    Belos::MinresSolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(belosList) );

    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Max number of CG iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = solver.solve();
    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<MT> actual_resids( numrhs );
    std::vector<MT> rhs_norm( numrhs );
    MV resid(map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -one, resid, one, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    }
    for ( int i=0; i<numrhs; i++) {
      MT actRes = actual_resids[i]/rhs_norm[i];
      if (proc_verbose) {
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      }
      if (actRes > tol) badRes = true;
    }

    success = (ret==Belos::Converged && !badRes);

    if (success) {
      if (proc_verbose)
        std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_bl_cg_hb.cpp
