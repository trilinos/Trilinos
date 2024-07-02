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
#include "BelosConfigDefs.hpp"
#include "BelosLinearMultiShiftProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosEpetraUtils.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[]) {
  //
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  //
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool verbose = false;
  bool success = true;
  try {
    bool proc_verbose = false;
    bool use_shift = true;
    int frequency = -1;    // how often residuals are printed by solver
    int init_numrhs = 3;   // how many right-hand sides get solved first
    int aug_numrhs = 2;   // how many right-hand sides are augmented to the first group
    int maxrestarts = 15;  // number of restarts allowed
    int length = 100;
    int init_blocksize = 2;// blocksize used for the initial pseudo-block GMRES solve
    int aug_blocksize = 1; // blocksize used for the augmented pseudo-block GMRES solve
    int maxiters = -1;     // maximum iterations allowed
    std::string filename("orsirr1.hb");
    MT tol = 1.0e-5;       // relative residual tolerance
    MT aug_tol = 1.0e-5;   // relative residual tolerance for augmented system

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("use-shift", "no-shift", &use_shift, "Whether shift should be used.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("aug-tol",&aug_tol,"Relative residual tolerance used by GMRES solver for augmented systems.");
    cmdp.setOption("init-num-rhs",&init_numrhs,"Number of right-hand sides to be initially solved for.");
    cmdp.setOption("aug-num-rhs",&aug_numrhs,"Number of right-hand sides augmenting the initial solve.");
    cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
    cmdp.setOption("block-size",&init_blocksize,"Block size used by GMRES for the initial solve.");
    cmdp.setOption("aug-block-size",&aug_blocksize,"Block size used by GMRES for the augmented solve.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    cmdp.setOption("subspace-size",&length,"Dimension of Krylov subspace used by GMRES.");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose
    //
    // Get the problem
    //
    int MyPID;
    RCP<Epetra_CrsMatrix> A;
    int return_val =Belos::Util::createEpetraProblem(filename,NULL,&A,NULL,NULL,&MyPID);
    const Epetra_Map &Map = A->RowMap();
    if(return_val != 0) return return_val;
    proc_verbose = verbose && (MyPID==0); /* Only print on zero processor */
    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = Map.NumGlobalElements();
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Num Blocks", length );                 // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", init_blocksize );         // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Deflation Quorum", init_blocksize  );  // Number of converged linear systems before deflation
    belosList.set( "Timer Label", "Belos Init" );          // Label used by the timers in this solver
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    //
    // *****Construct solution std::vector and random right-hand-sides *****
    //
    RCP<Epetra_MultiVector> initX = rcp( new Epetra_MultiVector(Map, init_numrhs) );
    RCP<Epetra_MultiVector> initB = rcp( new Epetra_MultiVector(Map, init_numrhs) );

    Belos::LinearMultiShiftProblem<double,MV,OP> initProblem( A, initX, initB );
    initProblem.setLabel("Belos Init");
    initProblem.setShift( use_shift );
    bool set = initProblem.setProblem(); 

    initX->Random();
    initProblem.applyOp( *initX, *initB );
    initX->PutScalar( 0.0 );
    set = initProblem.setProblem( initX, initB );

    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Initial Belos::LinearMultiShiftProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // *******************************************************************
    // *********************Perform initial solve*************************
    // *******************************************************************
    //
    Teuchos::RCP< Belos::SolverManager<double,MV,OP> > initSolver
      = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>( rcp(&initProblem,false), rcp(&belosList,false) ) );
    //
    // Perform solve
    //
    Belos::ReturnType ret = initSolver->solve();

    //
    // Compute actual residuals.
    //
    std::vector<int> init_index( init_numrhs );
    for (int i=0; i<init_numrhs; ++i)
      init_index[i] = i;
    initProblem.setLSIndex( init_index );

    bool badRes = false;
    std::vector<double> actual_resids( init_numrhs );
    std::vector<double> rhs_norm( init_numrhs );
    Epetra_MultiVector initR( Map, init_numrhs );
    initProblem.applyOp( *initX, initR );
    MVT::MvAddMv( -1.0, initR, 1.0, *initB, initR );
    MVT::MvNorm( initR, actual_resids );
    MVT::MvNorm( *initB, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for (int i=0; i<init_numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    if (ret != Belos::Converged || badRes==true) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Initial solve did not converge to solution!" << std::endl;
      return -1;
    }

    //
    // ***************Construct augmented linear system****************
    //
    RCP<Epetra_MultiVector> augX = rcp( new Epetra_MultiVector(Map, init_numrhs+aug_numrhs) );
    RCP<Epetra_MultiVector> augB = rcp( new Epetra_MultiVector(Map, init_numrhs+aug_numrhs) );
    Belos::LinearMultiShiftProblem<double,MV,OP> augProblem( A, augX, augB );
    augProblem.setLabel("Belos Aug");
    augProblem.setShift( use_shift );
    set = augProblem.setProblem();

    if (aug_numrhs) {
      augX->Random();
      augProblem.applyOp( *augX, *augB );
      augX->PutScalar( 0.0 );
    }

    // Copy previous linear system into
    RCP<Epetra_MultiVector> tmpX = rcp( new Epetra_MultiVector( Epetra_DataAccess::View, *augX, 0, init_numrhs ) );
    RCP<Epetra_MultiVector> tmpB = rcp( new Epetra_MultiVector( Epetra_DataAccess::View, *augB, 0, init_numrhs ) );
    tmpX->Scale( 1.0, *initX );
    tmpB->Scale( 1.0, *initB );

    set = augProblem.setProblem( augX, augB );
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Augmented Belos::LinearMultiShiftProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // *******************************************************************
    // *******************Perform augmented solve*************************
    // *******************************************************************
    //
    belosList.set( "Block Size", aug_blocksize );                // Blocksize to be used by iterative solver
    belosList.set( "Convergence Tolerance", aug_tol );           // Relative convergence tolerance requested
    belosList.set( "Deflation Quorum", aug_blocksize );          // Number of converged linear systems before deflation
    belosList.set( "Timer Label", "Belos Aug" );          // Label used by the timers in this solver
    belosList.set( "Implicit Residual Scaling", "Norm of RHS" ); // Implicit residual scaling for convergence
    belosList.set( "Explicit Residual Scaling", "Norm of RHS" ); // Explicit residual scaling for convergence
    Teuchos::RCP< Belos::SolverManager<double,MV,OP> > augSolver
      = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>( rcp(&augProblem,false), rcp(&belosList,false) ) );
    //
    // Perform solve
    //
    ret = augSolver->solve();

    if (ret != Belos::Converged) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR: Augmented solver did not converge to solution!" << std::endl;
      return -1;
    }
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of initial right-hand sides: " << init_numrhs << std::endl;
      std::cout << "Number of augmented right-hand sides: " << aug_numrhs << std::endl;
      std::cout << "Number of restarts allowed: " << maxrestarts << std::endl;
      std::cout << "Length of Krylov subspace: " << length <<std::endl;
      std::cout << "Max number of Gmres iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      if (aug_tol != tol)
        std::cout << "Relative residual tolerance for augmented systems: " << aug_tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Compute actual residuals.
    //
    int total_numrhs = init_numrhs + aug_numrhs;
    std::vector<int> total_index( total_numrhs );
    for (int i=0; i<total_numrhs; ++i)
      total_index[i] = i;
    augProblem.setLSIndex( total_index );

    badRes = false;
    actual_resids.resize( total_numrhs );
    rhs_norm.resize( total_numrhs );
    Epetra_MultiVector augR( Map, total_numrhs );
    augProblem.apply( *augX, augR );
    MVT::MvAddMv( -1.0, augR, 1.0, *augB, augR );
    MVT::MvNorm( augR, actual_resids );
    MVT::MvNorm( *augB, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<total_numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol ) badRes = true;
      }
    }
    if (ret!=Belos::Converged || badRes==true) {
      success = false;
      if (proc_verbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    } else {
      success = true;
      if (proc_verbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
} // end test_pseudo_gmres_hb.cpp
