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
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosFixedPointSolMgr.hpp"
#include "BelosEpetraUtils.h"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[]) {
  //
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  //
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef double                          ST;
  typedef Epetra_Operator                 OP;
  typedef Epetra_MultiVector              MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  
  bool success = false;
  bool verbose = false;
  try {
    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;  // total number of right-hand sides to solve for
    int maxiters = -1;  // maximum number of iterations for solver to use
    std::string filename("bcsstk14.hb");
    double tol = 1.0e-5;  // relative residual tolerance
    bool precond = false; // use diagonal preconditioner
    bool leftPrecond = false; // if preconditioner is used, left or right?


    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used byfixed point solver.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem).");
    cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
    cmdp.setOption("left","right",&leftPrecond,"Use a left/right preconditioner.");
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
    RCP<Epetra_MultiVector> B, X;
    int return_val =Belos::Util::createEpetraProblem(filename,NULL,&A,&B,&X,&MyPID);
    if(return_val != 0) return return_val;
    proc_verbose = ( verbose && (MyPID==0) );

    // Scale down the rhs
    std::vector<double> norm( MVT::GetNumberVecs(*B));
    MVT::MvNorm(*B,norm);
    for(int i=0; i< MVT::GetNumberVecs(*B); i++)
      norm[i] = 1.0 / norm[i];
    MVT::MvScale(*B,norm);

    // Drastically bump the diagonal and then rescale so FP can converge
    Epetra_Vector diag(A->RowMap());
    A->ExtractDiagonalCopy(diag);
    for(int i=0; i<diag.MyLength(); i++)
      diag[i]=diag[i]*1e4;
    A->ReplaceDiagonalValues(diag);

    for(int i=0; i<diag.MyLength(); i++)
      diag[i] = 1.0/sqrt(diag[i]);
    A->LeftScale(diag);
    A->RightScale(diag);

    //
    // Solve using Belos
    //
  
    //
    // *****Construct initial guess and random right-hand-sides *****
    //
    if (numrhs != 1) {
      X = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
      MVT::MvRandom( *X );
      B = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
      OPT::Apply( *A, *X, *B );
      MVT::MvInit( *X, 0.0 );
    }
    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->GlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    //
    // Construct an linear problem instance.
    //
    Belos::LinearProblem<double,MV,OP> problem( A, X, B );
    // diagonal preconditioner
    if (precond) {
      Epetra_Vector diagonal(A->RowMap());
      A->ExtractDiagonalCopy(diagonal);

      int NumMyElements    = diagonal.Map().NumMyElements();
      Teuchos::ArrayView<const int> MyGlobalElements = Teuchos::ArrayView< const int >(diagonal.Map().MyGlobalElements(), NumMyElements);
      RCP<Epetra_CrsMatrix> invDiagMatrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RowMap(), 1, true));

      for (Teuchos_Ordinal i=0; i<NumMyElements; ++i) {
        diagonal[i] = 1.0 / diagonal[i];
        invDiagMatrix->InsertGlobalValues(MyGlobalElements[i],
                                          1,
                                          &diagonal[i],
                                          &MyGlobalElements[i] );
      }
      invDiagMatrix->FillComplete();

      if (leftPrecond)
        problem.setLeftPrec(invDiagMatrix);
      else
        problem.setRightPrec(invDiagMatrix);
    }
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // Create an iterative solver manager.
    //
    RCP< Belos::SolverManager<double,MV,OP> > newSolver
      = rcp( new Belos::FixedPointSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = newSolver->solve();
    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    Epetra_MultiVector resid(A->Map(), numrhs);
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
        std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
