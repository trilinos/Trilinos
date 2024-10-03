// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// Multiple right-hand-sides are created randomly.
// The initial guesses are all set to zero.

/* Originally convert test here: belos/epetra/test/BlockCG/test_pseudo_cg_hb.cpp */

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"

#include <Trilinos_Util.h>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MatrixIO.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

template <typename ScalarType>
int run(int argc, char *argv[])
{
  using ST = ScalarType;
  using OP = typename Tpetra::Operator<ST>;
  using MV = typename Tpetra::MultiVector<ST>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST>;
  
  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT = typename SCT::magnitudeType;

  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  
  using Teuchos::ParameterList;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession session(&argc, &argv, &std::cout);
  RCP<const Comm<int>> comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);
  
  bool success = false;
  bool verbose = false;

  try
  {
    // Get test parameters from command-line processor
    bool proc_verbose = false;
    bool leftprec = true;      // left preconditioning or right.
    int frequency = -1; // how often residuals are printed by solver
    int numrhs = 1;  // total number of right-hand sides to solve for
    int maxiters = -1; // maximum number of iterations for the solver to use
    std::string filename("bcsstk14.hb");
    MT tol = 1.0e-5;  // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem size).");

    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose) {
      frequency = -1;  // Reset frequency if verbosity is off
    }

    proc_verbose = ( verbose && (MyPID==0) );
    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    // Get the problem
    RCP<tcrsmatrix_t> A;
    Tpetra::Utils::readHBMatrix(filename,comm,A);
    RCP<const Tpetra::Map<> > map = A->getDomainMap();

    // Construct initial guess and random right-hand-sides
    RCP<MV> B, X;
    X = rcp(new MV(map, numrhs));
    MVT::MvRandom(*X);
    B = rcp(new MV(map, numrhs));
    OPT::Apply(*A, *X, *B);
    MVT::MvInit(*X, 0.0);

    // Create parameter list for the pseudo block CG solver manager
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements; // maximum number of iterations to run
    
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Estimate Condition Number", true );
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    
    // Construct a linear problem
    Belos::LinearProblem<ST,MV,OP> problem(A, X, B);
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // Create an iterative solver manager.
    Belos::PseudoBlockCGSolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(belosList));

    // Start the block CG iteration
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of pseudo block CG iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    
    // Perform solve
    Belos::ReturnType ret = solver.solve();
    
    // Compute actual residuals.
    bool badRes = false;
    std::vector<MT> actual_resids( numrhs );
    std::vector<MT> rhs_norm( numrhs );
    MV resid(map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout << "---------- Actual Residuals (normalized) ----------" << std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        ST actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
      std::cout << std::endl << "Condition Estimate: " << solver.getConditionEstimate() << std::endl;
    }

    success = ret==Belos::Converged && !badRes;

    if (success) {
      if (proc_verbose)
        std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_pseudo_cg_hb.cpp

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);

  // wrapped with a check: CMake option Trilinos_ENABLE_FLOAT=ON
  // return run<float>(argc,argv);
}

