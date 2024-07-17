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
// The right-hand-side from the problem is being used instead of multiple
// random right-hand-sides.  The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//

// Tpetra
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"

template <typename ScalarType>
int run(int argc, char *argv[]) {
  //
  Teuchos::GlobalMPISession session(&argc, &argv, nullptr);
  //
  using ST = typename Tpetra::Vector<ScalarType>::scalar_type;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;

  using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT  = typename SCT::magnitudeType;
  using MGT = typename Teuchos::ScalarTraits<MT>;

  using tmap_t       = Tpetra::Map<LO,GO,NT>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  bool success = false;
  bool verbose = false;

  try {

    const auto comm = Tpetra::getDefaultComm();
    const int myPID = comm->getRank();

    bool badRes = false;
    bool procVerbose = false;
    bool pseudo = false;   // use pseudo block Gmres to solve this linear system.
    int frequency = -1;
    int blockSize = 1;
    int numrhs = 1;
    int maxRestarts = 15; // number of restarts allowed
    int maxIters = -1;    // maximum number of iterations allowed per linear system
    std::string filename("orsirr1.hb");
    std::string ortho("DGKS");
    MT tol = 1.0e-5;  // relative residual tolerance
    MT compTol = 10*MGT::prec();

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("pseudo","regular",&pseudo,"Use pseudo-block Gmres to solve the linear systems.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("ortho",&ortho,"Orthogonalization routine used by Gmres solver.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by Gmres solver.");
    cmdp.setOption("max-restarts",&maxRestarts,"Maximum number of restarts allowed for Gmres solver.");
    cmdp.setOption("blockSize",&blockSize,"Block size used by Gmres.");
    cmdp.setOption("maxIters",&maxIters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    procVerbose = ( verbose && (myPID==0) );

    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( comm );
    RCP<tcrsmatrix_t> A = reader.readFromFile( filename );
    RCP<const tmap_t> Map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp( new MV(Map,numrhs) );
    MVT::MvRandom( *X );
    B = rcp( new MV(Map,numrhs) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );

    // ********Other information used by block solver***********
    // *****************(can be user specified)******************

    const int NumGlobalElements = B->getGlobalLength();
    if (maxIters == -1)
      maxIters = NumGlobalElements/blockSize - 1; // maximum number of iterations to run

    ParameterList belosList;
    belosList.set( "Num Blocks", maxIters );               // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", blockSize );              // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxIters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxRestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Orthogonalization", ortho );           // Orthogonalization routine
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::StatusTestDetails );
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

    // *******************************************************************
    // *************Start the block Gmres iteration*************************
    // *******************************************************************

    RCP< Belos::SolverManager<ST,MV,OP> > solver;
    if (pseudo)
      solver = rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>( rcp(&problem,false), rcp(&belosList,false) ) );
    else
      solver = rcp( new Belos::BlockGmresSolMgr<ST,MV,OP>( rcp(&problem,false), rcp(&belosList,false) ) );

    // Perform solve
    Belos::ReturnType ret = solver->solve();

    if (ret!=Belos::Converged) {
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
      return -1;
    }

    // Get the number of iterations for this solve.
    int numIters = solver->getNumIters();
    if (procVerbose)
      std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;

    // Compute actual residuals.
    std::vector<MT> actual_resids( numrhs );
    std::vector<MT> rhs_norm( numrhs );
    MV resid(Map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
      }
    }

    // -----------------------------------------------------------------
    // Resolve the first problem by just resetting the solver manager.
    // -----------------------------------------------------------------

    X->putScalar( 0.0 );
    solver->reset( Belos::Problem );

    // Perform solve (again)
    ret = solver->solve();

    // Get the number of iterations for this solve.
    numIters = solver->getNumIters();
    if (procVerbose)
      std::cout << "Number of iterations performed for this solve (manager reset): " << numIters << std::endl;

    if (ret!=Belos::Converged) {
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
      return -1;
    }

    // Compute actual residuals.
    std::vector<MT> actual_resids2( numrhs );
    MV resid2(Map, numrhs);
    OPT::Apply( *A, *X, resid2 );
    MVT::MvAddMv( -1.0, resid2, 1.0, *B, resid2 );
    MVT::MvNorm( resid2, actual_resids );
    MVT::MvAddMv( -1.0, resid2, 1.0, resid, resid2 );
    MVT::MvNorm( resid2, actual_resids2 );
    MVT::MvNorm( *B, rhs_norm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (manager reset) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
        if ( actual_resids2[i]/rhs_norm[i] > compTol ) {
          badRes = true;
          std::cout << "Resolve residual vector is too different from first solve residual vector: " << actual_resids2[i]/rhs_norm[i] << std::endl;
        }
      }
    }

    // -----------------------------------------------------------------
    // Resolve the first problem by resetting the linear problem.
    // -----------------------------------------------------------------

    X->putScalar( 0.0 );
    set = problem.setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    solver->setProblem( rcp(&problem,false) );

    // Perform solve (again)
    ret = solver->solve();

    // Get the number of iterations for this solve.
    numIters = solver->getNumIters();
    if (procVerbose)
      std::cout << "Number of iterations performed for this solve (manager setProblem()): " << numIters << std::endl;

    if (ret!=Belos::Converged) {
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
      return -1;
    }

    // Compute actual residuals.
    OPT::Apply( *A, *X, resid2 );
    MVT::MvAddMv( -1.0, resid2, 1.0, *B, resid2 );
    MVT::MvNorm( resid2, actual_resids );
    MVT::MvAddMv( -1.0, resid2, 1.0, resid, resid2 );
    MVT::MvNorm( resid2, actual_resids2 );
    MVT::MvNorm( *B, rhs_norm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (manager setProblem()) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
        if ( actual_resids2[i]/rhs_norm[i] > compTol ) {
          badRes = true;
          std::cout << "Resolve residual vector is too different from first solve residual vector: " << actual_resids2[i]/rhs_norm[i] << std::endl;
        }
      }
    }

    // ----------------------------------------------------------------------------------
    // Resolve the first problem by resetting the solver manager and changing the labels.
    // ----------------------------------------------------------------------------------

    ParameterList belosList2;
    belosList2.set( "Timer Label", "Belos Resolve w/ New Label" );   // Change timer labels.
    solver->setParameters( rcp( &belosList2, false ) );

    problem.setLabel( "Belos Resolve w/ New Label" );
    X->putScalar( 0.0 );
    solver->reset( Belos::Problem );

    // Perform solve (again)
    ret = solver->solve();

    // Get the number of iterations for this solve.
    numIters = solver->getNumIters();
    if (procVerbose)
      std::cout << "Number of iterations performed for this solve (label reset): " << numIters << std::endl;

    if (ret!=Belos::Converged) {
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
      return -1;
    }

    // Compute actual residuals.
    OPT::Apply( *A, *X, resid2 );
    MVT::MvAddMv( -1.0, resid2, 1.0, *B, resid2 );
    MVT::MvNorm( resid2, actual_resids );
    MVT::MvAddMv( -1.0, resid2, 1.0, resid, resid2 );
    MVT::MvNorm( resid2, actual_resids2 );
    MVT::MvNorm( *B, rhs_norm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (label reset) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
        if ( actual_resids2[i]/rhs_norm[i] > compTol ) {
          badRes = true;
          std::cout << "Resolve residual vector is too different from first solve residual vector: " << actual_resids2[i]/rhs_norm[i] << std::endl;
        }
      }
    }

    // -------------------------------------------------------------
    // Construct a second unpreconditioned linear problem instance.
    // -------------------------------------------------------------

    RCP<MV> X2 = MVT::Clone(*X, numrhs);
    MVT::MvInit( *X2, 0.0 );
    Belos::LinearProblem<ST,MV,OP> problem2( A, X2, B );
    problem2.setLabel("Belos Resolve");
    set = problem2.setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // *******************************************************************
    // *************Start the block Gmres iteration*************************
    // *******************************************************************

    // Create the solver without either the problem or parameter list.
    if (pseudo)
      solver = rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>() );
    else
      solver = rcp( new Belos::BlockGmresSolMgr<ST,MV,OP>() );

    // Set the problem after the solver construction.
    solver->setProblem( rcp( &problem2, false ) );

    // Get the valid list of parameters from the solver and print it.
    RCP<const ParameterList> validList = solver->getValidParameters();
    if (procVerbose) {
      if (pseudo)
        std::cout << std::endl << "Valid parameters from the pseudo-block Gmres solver manager:" << std::endl;
      else
        std::cout << std::endl << "Valid parameters from the block Gmres solver manager:" << std::endl;

      std::cout << *validList << std::endl;
    }

    // Set the parameter list after the solver construction.
    belosList.set( "Timer Label", "Belos Resolve" );         // Set timer label to discern between the two solvers.
    solver->setParameters( rcp( &belosList, false ) );

    // Perform solve
    ret = solver->solve();

    // Get the number of iterations for this solve.
    numIters = solver->getNumIters();
    if (procVerbose)
      std::cout << "Number of iterations performed for this solve (new solver): " << numIters << std::endl;

    // Compute actual residuals.
    OPT::Apply( *A, *X2, resid2 );
    MVT::MvAddMv( -1.0, resid2, 1.0, *B, resid2 );
    MVT::MvNorm( resid2, actual_resids );
    MVT::MvAddMv( -1.0, resid2, 1.0, resid, resid2 );
    MVT::MvNorm( resid2, actual_resids2 );
    MVT::MvNorm( *B, rhs_norm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (new solver) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
        if ( actual_resids2[i]/rhs_norm[i] > compTol ) {
          badRes = true;
          std::cout << "Resolve residual vector is too different from first solve residual vector: " << actual_resids2[i]/rhs_norm[i] << std::endl;
        }
      }
    }

    // **********Print out information about problem*******************
    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blockSize << std::endl;
      std::cout << "Max number of Gmres iterations: " << maxIters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }

    success = ret==Belos::Converged && !badRes;

    if (success) {
      if (procVerbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    } else {
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // run

int main(int argc, char *argv[]) {
  // run with different ST
  return run<double>(argc,argv);
  // run<float>(argc,argv); // FAILS
}
