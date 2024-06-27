// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// NOTE: No preconditioner is used in this example.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"

#include <Teuchos_GlobalMPISession.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

template<typename ScalarType>
int run(int argc, char *argv[]) {
  using Teuchos::CommandLineProcessor;
  using Teuchos::GlobalMPISession;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::tuple;

  using ST  = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO  = typename Tpetra::Vector<>::local_ordinal_type;
  using GO  = typename Tpetra::Vector<>::global_ordinal_type;
  using NT  = typename Tpetra::Vector<>::node_type;
  
  using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MAP = typename Tpetra::Map<LO,GO,NT>;

  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MAT = Tpetra::CrsMatrix<ST,LO,GO,NT>;
  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT  = typename SCT::magnitudeType;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  bool success = false;
  bool verbose = false;

  try
  {
    // Get test parameters from command-line processor
    //
    bool debug = false, proc_verbose = false;
    int frequency = -1;        // frequency of status test output.
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxiters = -1;         // maximum number of iterations allowed per linear system
    MT tol = 1.0e-10;           // relative residual tolerance

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information from the solver.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose
    // **********************************************************************
    // ******************Set up the problem to be solved*********************
    // construct diagonal matrix
    const GO numGlobalElements = 100;
    const int m = 4; // number of negative eigenvalues
    // Create diagonal matrix with n-m positive and m negative eigenvalues.
    RCP<const MAP> map = rcp(new MAP(numGlobalElements, 0, comm));

    RCP<MAT> A = rcp( new MAT(map, 1) );

    for ( size_t k = 0; k < map->getLocalNumElements(); k++ )
    {
      GO gid = map->getGlobalElement(k);
      const ST val = static_cast<ST> (2 * (gid-m) + 1);
      A->insertGlobalValues(gid, tuple<GO> (gid), tuple<ST> (val));
    }

    A->fillComplete();
    TEUCHOS_ASSERT(A->isStorageOptimized());

    // create initial guess and right-hand side
    RCP<MV> vecX = rcp( new MV( map, numrhs ) );

    RCP<MV> vecB = rcp( new MV( map, numrhs ) );
    // **********************************************************************
    
    proc_verbose = ( verbose && (comm->getRank()==0) );  /* Only print on the zero processor */

    RCP<MV> X;
    RCP<MV> B;
    // Check to see if the number of right-hand sides is the same as requested.
    if (numrhs>1) {
      X = rcp( new MV( map, numrhs ) );
      B = rcp( new MV( map, numrhs ) );
      X->randomize();
      OPT::Apply( *A, *X, *B );
      X->putScalar( 0.0 );
    }
    else {
      X = rcp_implicit_cast<MV>(vecX);
      B = rcp_implicit_cast<MV>(vecB);
      B->putScalar( 1.0 );
    }
    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    if (maxiters == -1)
      maxiters = numGlobalElements - 1; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );        // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );          // Relative convergence tolerance requested
    belosList.set( "Assert Positive Definiteness", false ); // Explicitly don't enforce positive definiteness

    int verbosity = Belos::Errors + Belos::Warnings;
    if (verbose) {
      verbosity += Belos::TimingDetails + Belos::StatusTestDetails;
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    if (debug) {
      verbosity += Belos::Debug;
    }
    belosList.set( "Verbosity", verbosity );
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
    // ****************Start the CG iteration*************************
    // *******************************************************************
    //
    // Create an iterative solver manager.
    RCP<Belos::SolverManager<ST,MV,OP> > newSolver
      = rcp( new Belos::PseudoBlockCGSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << numGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = newSolver->solve();
    //
    // Get the number of iterations for this solve.
    //
    int numIters = newSolver->getNumIters();
    if (proc_verbose)
      std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;
    //
    // Compute actual residuals.
    //
    
    bool badRes = false;

    std::vector<ST> actual_resids( numrhs );
    std::vector<ST> rhs_norm( numrhs );
    MV resid(map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        ST actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
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

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);
} // end test_pseudo_cg_indefinite.cpp
