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
//
// This test is for testing the deflation in the pseudo-block Gmres solver.
// One set of linear systems is solved and then augmented with additional
// linear systems and resolved.  The already solved linear systems should be
// deflated immediately, leaving only the augmented systems to be solved.
//
//

// Tpetra
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosTpetraTestFramework.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosLinearMultiShiftProblem.hpp"

template <typename ScalarType>
int run(int argc, char *argv[]) {

  using ST = typename Tpetra::CrsMatrix<ScalarType>::scalar_type;
  using LO = typename Tpetra::CrsMatrix<>::local_ordinal_type;
  using GO = typename Tpetra::CrsMatrix<>::global_ordinal_type;
  using NT = typename Tpetra::CrsMatrix<>::node_type;

  using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MT  = typename Teuchos::ScalarTraits<ST>::magnitudeType;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;

  using tmap_t       = Tpetra::Map<LO,GO,NT>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  Teuchos::GlobalMPISession session(&argc, &argv, nullptr);

  bool verbose = false;
  bool success = true;

  try {

    const auto comm = Tpetra::getDefaultComm();
    const int myPID = comm->getRank();

    bool procVerbose = false;
    bool useShift = true;
    int frequency = -1;    // how often residuals are printed by solver
    int initNumRHS = 3;   // how many right-hand sides get solved first
    int augNumRHS = 2;   // how many right-hand sides are augmented to the first group
    int maxRestarts = 15;  // number of restarts allowed
    int length = 100;
    int initBlockSize = 2;// blocksize used for the initial pseudo-block GMRES solve
    int augBlockSize = 1; // blocksize used for the augmented pseudo-block GMRES solve
    int maxIters = -1;     // maximum iterations allowed
    std::string filename("orsirr1.hb");
    MT tol = 1.0e-5;       // relative residual tolerance
    MT aug_tol = 1.0e-5;   // relative residual tolerance for augmented system

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("use-shift", "no-shift", &useShift, "Whether shift should be used.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("aug-tol",&aug_tol,"Relative residual tolerance used by GMRES solver for augmented systems.");
    cmdp.setOption("init-num-rhs",&initNumRHS,"Number of right-hand sides to be initially solved for.");
    cmdp.setOption("aug-num-rhs",&augNumRHS,"Number of right-hand sides augmenting the initial solve.");
    cmdp.setOption("max-restarts",&maxRestarts,"Maximum number of restarts allowed for GMRES solver.");
    cmdp.setOption("block-size",&initBlockSize,"Block size used by GMRES for the initial solve.");
    cmdp.setOption("aug-block-size",&augBlockSize,"Block size used by GMRES for the augmented solve.");
    cmdp.setOption("max-iters",&maxIters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    cmdp.setOption("subspace-size",&length,"Dimension of Krylov subspace used by GMRES.");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }
    procVerbose = verbose && (myPID==0); /* Only print on zero processor */

    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( comm );
    RCP<tcrsmatrix_t> A = reader.readFromFile( filename );
    RCP<const tmap_t> Map = A->getDomainMap();

    // ********Other information used by block solver***********
    // *****************(can be user specified)******************

    const int numGlobalElements = Map->getGlobalNumElements();
    if (maxIters == -1)
      maxIters = numGlobalElements - 1; // maximum number of iterations to run

    ParameterList belosList;
    belosList.set( "Num Blocks", length );                 // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", initBlockSize );         // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxIters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxRestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Deflation Quorum", initBlockSize  );  // Number of converged linear systems before deflation
    belosList.set( "Timer Label", "Belos Init" );          // Label used by the timers in this solver
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

    // *****Construct solution std::vector and random right-hand-sides *****
    RCP<MV> initX = rcp( new MV(Map, initNumRHS) );
    RCP<MV> initB = rcp( new MV(Map, initNumRHS) );

    Belos::LinearMultiShiftProblem<ST,MV,OP> initProblem( A, initX, initB );
    initProblem.setLabel("Belos Init");
    initProblem.setShift( useShift );
    bool set = initProblem.setProblem();

    MVT::MvRandom( *initX );
    initProblem.applyOp( *initX, *initB );
    initX->putScalar( 0.0 );
    set = initProblem.setProblem( initX, initB );

    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Initial Belos::LinearMultiShiftProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // *******************************************************************
    // *********************Perform initial solve*************************
    // *******************************************************************

    RCP< Belos::SolverManager<ST,MV,OP> > initSolver
      = rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>( rcp(&initProblem,false), rcp(&belosList,false) ) );

    // Perform solve
    Belos::ReturnType ret = initSolver->solve();

    // Compute actual residuals.
    std::vector<int> initIndex( initNumRHS );
    for (int i=0; i<initNumRHS; ++i)
      initIndex[i] = i;
    initProblem.setLSIndex( initIndex );

    bool badRes = false;
    std::vector<ST> actualResids( initNumRHS );
    std::vector<ST> rhsNorm( initNumRHS );
    MV initR( Map, initNumRHS );
    initProblem.applyOp( *initX, initR );
    MVT::MvAddMv( -1.0, initR, 1.0, *initB, initR );
    MVT::MvNorm( initR, actualResids );
    MVT::MvNorm( *initB, rhsNorm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for (int i=0; i<initNumRHS; i++) {
        ST actRes = actualResids[i]/rhsNorm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    if (ret != Belos::Converged || badRes==true) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Initial solve did not converge to solution!" << std::endl;
      return -1;
    }

    // ***************Construct augmented linear system****************

    RCP<MV> augX = rcp( new MV(Map, initNumRHS+augNumRHS) );
    RCP<MV> augB = rcp( new MV(Map, initNumRHS+augNumRHS) );
    Belos::LinearMultiShiftProblem<ST,MV,OP> augProblem( A, augX, augB );
    augProblem.setLabel("Belos Aug");
    augProblem.setShift( useShift );
    set = augProblem.setProblem();

    if (augNumRHS) {
      MVT::MvRandom( *augX );
      augProblem.applyOp( *augX, *augB );
      augX->putScalar( 0.0 );
    }

    // Copy previous linear system into
    RCP<MV> tmpX = rcp( new MV( *augX ) );
    RCP<MV> tmpB = rcp( new MV( *augB ) );
    tmpX->scale( 1.0, *augX );
    tmpB->scale( 1.0, *augB );

    set = augProblem.setProblem( augX, augB );
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Augmented Belos::LinearMultiShiftProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // *******************************************************************
    // *******************Perform augmented solve*************************
    // *******************************************************************

    belosList.set( "Block Size", augBlockSize );                // Blocksize to be used by iterative solver
    belosList.set( "Convergence Tolerance", aug_tol );           // Relative convergence tolerance requested
    belosList.set( "Deflation Quorum", augBlockSize );          // Number of converged linear systems before deflation
    belosList.set( "Timer Label", "Belos Aug" );          // Label used by the timers in this solver
    belosList.set( "Implicit Residual Scaling", "Norm of RHS" ); // Implicit residual scaling for convergence
    belosList.set( "Explicit Residual Scaling", "Norm of RHS" ); // Explicit residual scaling for convergence
    RCP< Belos::SolverManager<ST,MV,OP> > augSolver
      = rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>( rcp(&augProblem,false), rcp(&belosList,false) ) );

    // Perform solve
    ret = augSolver->solve();

    if (ret != Belos::Converged) {
      if (procVerbose)
        std::cout << std::endl << "ERROR: Augmented solver did not converge to solution!" << std::endl;
      return -1;
    }

    // **********Print out information about problem*******************

    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << numGlobalElements << std::endl;
      std::cout << "Number of initial right-hand sides: " << initNumRHS << std::endl;
      std::cout << "Number of augmented right-hand sides: " << augNumRHS << std::endl;
      std::cout << "Number of restarts allowed: " << maxRestarts << std::endl;
      std::cout << "Length of Krylov subspace: " << length <<std::endl;
      std::cout << "Max number of Gmres iterations: " << maxIters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      if (aug_tol != tol) {
        std::cout << "Relative residual tolerance for augmented systems: " << aug_tol << std::endl;
      }
      std::cout << std::endl;
    }

    // Compute actual residuals.
    int totalNumRHS = initNumRHS + augNumRHS;
    std::vector<int> totalIndex( totalNumRHS );
    for (int i=0; i<totalNumRHS; ++i)
      totalIndex[i] = i;
    augProblem.setLSIndex( totalIndex );

    badRes = false;
    actualResids.resize( totalNumRHS );
    rhsNorm.resize( totalNumRHS );
    MV augR( Map, totalNumRHS );
    augProblem.apply( *augX, augR );
    MVT::MvAddMv( -1.0, augR, 1.0, *augB, augR );
    MVT::MvNorm( augR, actualResids );
    MVT::MvNorm( *augB, rhsNorm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<totalNumRHS; i++) {
        ST actRes = actualResids[i]/rhsNorm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol ) badRes = true;
      }
    }
    if (ret!=Belos::Converged || badRes==true) {
      success = false;
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    } else {
      success = true;
      if (procVerbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
} // run

int main(int argc, char *argv[]) {
  // run with different ST
  return run<double>(argc,argv);
  // run<float>(argc,argv); // FAILS
}
