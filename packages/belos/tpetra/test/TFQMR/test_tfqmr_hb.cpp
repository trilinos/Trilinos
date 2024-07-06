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
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>

// Teuchos
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTFQMRSolMgr.hpp"
#include "BelosPseudoBlockTFQMRSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"

template<typename ScalarType>
int run(int argc, char *argv[]) {  
  // Teuchos
  using SCT = typename Teuchos::ScalarTraits<ScalarType>;
  using MT = typename SCT::magnitudeType;
  
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Tpetra
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<>::node_type;

  using OP = Tpetra::Operator<ST,LO,GO,NT>;
  using MV = Tpetra::MultiVector<ST,LO,GO,NT>;

  using tcrsmatrix_t = Tpetra::CrsMatrix<ST>;

  // Belos
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;

  // MPI session
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);
  
  const ST one  = SCT::one();
  bool success = false;
  bool verbose = false;
  
  try
  {
    // Get test parameters from command-line processor
    bool procVerbose = false;
    bool explicitTest = false;
    bool compRecursive = false;
    bool pseudo = false;
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;  // total number of right-hand sides to solve for
    int maxiters = -1;  // maximum number of iterations for solver to use
    std::string filename("orsirr1.hb");
    MT tol = 1.0e-5; // relative residual tolerance
    
    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("explicit","implicit-only",&explicitTest,"Compute explicit residuals.");
    cmdp.setOption("recursive","native",&compRecursive,"Compute recursive residuals.");
    cmdp.setOption("pseudo","not-pseudo",&pseudo,"Use pseudo-block TFQMR solver.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by TFQMR or pseudo-block TFQMR solver.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( comm );
    RCP<tcrsmatrix_t> A = reader.readFromFile( filename );
    RCP<const Tpetra::Map<> > map = A->getDomainMap();
    
    // Create initial vectors
    RCP<MV> B, X;
    X = rcp( new MV(map,numrhs) );
    MVT::MvRandom( *X );
    B = rcp( new MV(map,numrhs) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );
    procVerbose = ( verbose && (MyPID==0) );
    
    // Solve using Belos
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements; // maximum number of iterations to run
    
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (explicitTest)
    {
      belosList.set( "Explicit Residual Test", true );       // Scale by norm of right-hand side vector."
      belosList.set( "Explicit Residual Scaling", "Norm of RHS" ); // Scale by norm of right-hand side vector."
    }
    if (compRecursive)
    {
      belosList.set( "Compute Recursive Residuals", true );
    }
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    
    // Construct an unpreconditioned linear problem instance.
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    bool set = problem.setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // Create an iterative solver manager.
    RCP< Belos::SolverManager<ST,MV,OP> > solver;
    if (pseudo)
      solver = rcp( new Belos::PseudoBlockTFQMRSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&belosList,false)) );
    else
      solver = rcp( new Belos::TFQMRSolMgr<ST,MV,OP>(rcp(&problem, false), rcp(&belosList, false)));
    
    // **********Print out information about problem*******************
    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of TFQMR iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }

    // Perform solve
    Belos::ReturnType ret = solver->solve();

    // Compute actual residuals.
    bool badRes = false;
    std::vector<MT> actualResids( numrhs );
    std::vector<MT> rhsNorm( numrhs );
    MV resid(map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -one, resid, one, *B, resid );
    MVT::MvNorm( resid, actualResids );
    MVT::MvNorm( *B, rhsNorm );
    
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        MT actRes = actualResids[i]/rhsNorm[i];
        if (procVerbose) {
          std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        }
        if (actRes > tol) badRes = true;
      }
    }
    
    success = ret==Belos::Converged && !badRes;

    if (success) {
      if (procVerbose)
        std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
    } else {
      if (procVerbose)
        std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);

  // wrapped with a check: CMake option Trilinos_ENABLE_FLOAT=ON
  // run<float>(argc,argv);
}
