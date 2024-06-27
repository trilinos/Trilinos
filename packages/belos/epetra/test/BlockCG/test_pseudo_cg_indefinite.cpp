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
#include "BelosEpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"

#include "Epetra_Map.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[]) {
  //
  int MyPID = 0;
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  MyPID = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
#endif
  //
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool success = false;
  bool verbose = false;

  try {
    bool debug = false, proc_verbose = false;
    int frequency = -1;        // frequency of status test output.
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxiters = -1;         // maximum number of iterations allowed per linear system
    MT tol = 1.0e-10;           // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information from the solver.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose
    // **********************************************************************
    // ******************Set up the problem to be solved*********************
    // construct diagonal matrix
    const int NumGlobalElements = 100;
    const int m = 4; // number of negative eigenvalues

    // Create diagonal matrix with n-m positive and m negative eigenvalues.
    Epetra_Map epetraMap( NumGlobalElements, 0, Comm );

    Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix( Epetra_DataAccess::Copy, epetraMap, 1 ) );
    for ( int k=0; k<epetraMap.NumMyElements(); k++ )
    {
      int GIDk = epetraMap.GID(k);
      double val = 2*(GIDk-m) + 1;
      TEUCHOS_ASSERT_EQUALITY( 0, A->InsertGlobalValues( GIDk, 1, &val, &GIDk ) );
    }
    TEUCHOS_ASSERT_EQUALITY( 0, A->FillComplete() );
    TEUCHOS_ASSERT_EQUALITY( 0, A->OptimizeStorage() );

    // create initial guess and right-hand side
    Teuchos::RCP<Epetra_MultiVector> vecX = Teuchos::rcp( new Epetra_MultiVector( epetraMap, numrhs ) );

    Teuchos::RCP<Epetra_MultiVector> vecB = Teuchos::rcp( new Epetra_MultiVector( epetraMap, numrhs ) );
    // **********************************************************************
    proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */


    Teuchos::RCP<Epetra_MultiVector> X;
    Teuchos::RCP<Epetra_MultiVector> B;
    // Check to see if the number of right-hand sides is the same as requested.
    if (numrhs>1) {
      X = rcp( new Epetra_MultiVector( epetraMap, numrhs ) );
      B = rcp( new Epetra_MultiVector( epetraMap, numrhs ) );
      X->Random();
      OPT::Apply( *A, *X, *B );
      X->PutScalar( 0.0 );
    }
    else {
      X = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecX);
      B = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecB);
      B->PutScalar( 1.0 );
    }
    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
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
    Belos::LinearProblem<double,MV,OP> problem( A, X, B );
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
    RCP<Belos::SolverManager<double,MV,OP> > newSolver
      = rcp( new Belos::PseudoBlockCGSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
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
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    Epetra_MultiVector resid(epetraMap, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
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

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
