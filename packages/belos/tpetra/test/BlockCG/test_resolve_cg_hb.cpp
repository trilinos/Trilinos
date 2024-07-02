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

/* Originally convert test here: belos/epetra/test/BlockCG/test_resolve_cg_hb.cpp */

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"

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
  using MT  = typename SCT::magnitudeType;
  
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  
  using Teuchos::ParameterList;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Comm<int>> comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);

  bool success = false;
  bool verbose = false;
  
  try
  {
    // Get test parameters from command-line processor
    bool badRes = false;
    bool proc_verbose = false;
    bool pseudo = false;   // use pseudo block CG to solve this linear system.
    int frequency = -1;
    int blocksize = 1;
    int numrhs = 1;
    int maxiters = -1;    // maximum number of iterations allowed per linear system
    std::string filename("bcsstk14.hb");
    std::string ortho("DGKS");
    MT tol = 1.0e-5;  // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("pseudo","regular",&pseudo,"Use pseudo-block CG to solve the linear systems.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("ortho",&ortho,"Orthogonalization routine used by CG solver.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
    cmdp.setOption("blocksize",&blocksize,"Block size used by CG.");
    cmdp.setOption("maxiters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose){
      frequency = -1;  // reset frequency if test is not verbose
    }

    proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */
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
    
    // Other information used by block solver (can be user specified)
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
    
    //
    ParameterList belosList;
    belosList.set( "Num Blocks", maxiters );               // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
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
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // Start the block CG iteration
    RCP< Belos::SolverManager<ST,MV,OP> > solver;
    if (pseudo)
      solver = rcp( new Belos::PseudoBlockCGSolMgr<ST,MV,OP>( rcpFromRef(problem), rcpFromRef(belosList) ) );
    else
      solver = rcp( new Belos::BlockCGSolMgr<ST,MV,OP>( rcpFromRef(problem), rcpFromRef(belosList) ) );
    
    // Perform solve
    Belos::ReturnType ret = solver->solve();
    if (ret!=Belos::Converged) {
      if (proc_verbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
      return -1;
    }
    
    // Get the number of iterations for this solve.
    int numIters = solver->getNumIters();
    if (proc_verbose)
      std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;
    
    // Compute actual residuals.
    std::vector<MT> actual_resids( numrhs );
    std::vector<MT> rhs_norm( numrhs );
    MV resid(map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
      }
    }
    
    // Resolve the first problem by just resetting the solver manager.
    X->putScalar( 0.0 );
    solver->reset( Belos::Problem );
    
    // Perform solve (again)
    ret = solver->solve();
    
    // Get the number of iterations for this solve.
    numIters = solver->getNumIters();
    if (proc_verbose)
      std::cout << "Number of iterations performed for this solve (manager reset): " << numIters << std::endl;
    if (ret!=Belos::Converged) {
      if (proc_verbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
      return -1;
    }
    
    // Compute actual residuals.
    std::vector<MT> actual_resids2( numrhs );
    MV resid2(map, numrhs);
    OPT::Apply( *A, *X, resid2 );
    MVT::MvAddMv( -1.0, resid2, 1.0, *B, resid2 );
    MVT::MvNorm( resid2, actual_resids );
    MVT::MvAddMv( -1.0, resid2, 1.0, resid, resid2 );
    MVT::MvNorm( resid2, actual_resids2 );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (manager reset) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
        if ( actual_resids2[i] > SCT::prec() ) {
          badRes = true;
          std::cout << "Resolve residual vector is too different from first solve residual vector: " << actual_resids2[i] << std::endl;
        }
      }
    }
    
    // Resolve the first problem by resetting the solver manager and changing the labels.
    ParameterList belosList2;
    belosList2.set( "Timer Label", "Belos Resolve w/ New Label" );   // Change timer labels.
    solver->setParameters( Teuchos::rcp( &belosList2, false ) );

    problem.setLabel( "Belos Resolve w/ New Label" );
    X->putScalar( 0.0 );
    solver->reset( Belos::Problem );
    
    // Perform solve (again)
    ret = solver->solve();
    
    // Get the number of iterations for this solve.
    numIters = solver->getNumIters();
    if (proc_verbose)
      std::cout << "Number of iterations performed for this solve (label reset): " << numIters << std::endl;
    if (ret!=Belos::Converged) {
      if (proc_verbose)
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
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (label reset) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
        if ( actual_resids2[i] > SCT::prec() ) {
          badRes = true;
          std::cout << "Resolve residual vector is too different from first solve residual vector: " << actual_resids2[i] << std::endl;
        }
      }
    }

    // Construct a second unpreconditioned linear problem instance.
    RCP<MV> X2 = MVT::Clone(*X, numrhs);
    MVT::MvInit( *X2, 0.0 );
    Belos::LinearProblem<ST,MV,OP> problem2( A, X2, B );
    problem2.setLabel("Belos Resolve");
    set = problem2.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // Start the block CG iteration
    // Create the solver without either the problem or parameter list.
    if (pseudo)
      solver = rcp( new Belos::PseudoBlockCGSolMgr<ST,MV,OP>() );
    else
      solver = rcp( new Belos::BlockCGSolMgr<ST,MV,OP>() );
    
    // Set the problem after the solver construction.
    solver->setProblem( rcpFromRef(problem2) );

    // Get the valid list of parameters from the solver and print it.
    RCP<const ParameterList> validList = solver->getValidParameters();
    if (proc_verbose) {
      if (pseudo)
        std::cout << std::endl << "Valid parameters from the pseudo-block CG solver manager:" << std::endl;
      else
        std::cout << std::endl << "Valid parameters from the block CG solver manager:" << std::endl;

      std::cout << *validList << std::endl;
    }

    // Set the parameter list after the solver construction.
    belosList.set( "Timer Label", "Belos Resolve" );         // Set timer label to discern between the two solvers.
    solver->setParameters( rcp( &belosList, false ) );
    
    // Perform solve
    ret = solver->solve();
    
    // Get the number of iterations for this solve.
    numIters = solver->getNumIters();
    if (proc_verbose)
      std::cout << "Number of iterations performed for this solve (new solver): " << numIters << std::endl;
    
    // Compute actual residuals.
    OPT::Apply( *A, *X2, resid2 );
    MVT::MvAddMv( -1.0, resid2, 1.0, *B, resid2 );
    MVT::MvNorm( resid2, actual_resids );
    MVT::MvAddMv( -1.0, resid2, 1.0, resid, resid2 );
    MVT::MvNorm( resid2, actual_resids2 );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (new solver) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
        if ( actual_resids2[i] > SCT::prec() ) {
          badRes = true;
          std::cout << "Resolve residual vector is too different from first solve residual vector: " << actual_resids2[i] << std::endl;
        } 
      }
    }

    // Print out information about problem
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Max number of CG iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }

    success = ret==Belos::Converged && !badRes;

    if (success) {
      if (proc_verbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_resolve_cg_hb.cpp

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);

  // wrapped with a check: CMake option Trilinos_ENABLE_FLOAT=ON
  // return run<float>(argc,argv);
}

