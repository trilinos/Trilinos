// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a file, which can be in Harwell-Boeing (*.hb),
// Matrix Market (*.mtx), or triplet format (*.triU, *.triS).  The right-hand side
// from the problem, if it exists, will be used instead of multiple
// random right-hand-sides.  The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this example.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"
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

bool verbose = false;
bool success = true;
try {
bool proc_verbose = false;
  int frequency = -1;        // frequency of status test output.
  int blocksize = 1;         // blocksize
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  int maxsubspace = 50;      // maximum number of blocks the solver can use for the subspace
  int maxrestarts = 15;      // number of restarts allowed
  std::string filename("orsirr1.hb");
  MT tol = 1.0e-5;           // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blocksize,"Block size used by GMRES.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
  cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of blocks the solver can use for the subspace.");
  cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose
  //
  // Get the problem
  //
  RCP<Epetra_Map> Map;
  RCP<Epetra_CrsMatrix> A;
  RCP<Epetra_MultiVector> B, X;
  RCP<Epetra_Vector> vecB, vecX;
  EpetraExt::readEpetraLinearSystem(filename, Comm, &A, &Map, &vecX, &vecB);
  A->OptimizeStorage();
  proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

  // Check to see if the number of right-hand sides is the same as requested.
  if (numrhs>1) {
    X = rcp( new Epetra_MultiVector( *Map, numrhs ) );
    B = rcp( new Epetra_MultiVector( *Map, numrhs ) );
    X->Random();
    OPT::Apply( *A, *X, *B );
    X->PutScalar( 0.0 );
  }
  else {
    X = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecX);
    B = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecB);
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
  belosList.set( "Num Blocks", maxsubspace);             // Maximum number of blocks in Krylov factorization
  belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  if (verbose) {
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
		   Belos::TimingDetails + Belos::StatusTestDetails );
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  else
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
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
  // *************Start the block Gmres iteration*************************
  // *******************************************************************
  //
  // Create an iterative solver manager.
  RCP< Belos::SolverManager<double,MV,OP> > newSolver
    = rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)) );

  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Max number of restarts allowed: " << maxrestarts << std::endl;
    std::cout << "Max number of Gmres iterations per restart cycle: " << maxiters << std::endl;
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
  Epetra_MultiVector resid(*Map, numrhs);
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

if (ret!=Belos::Converged || badRes) {
  success = false;
  if (proc_verbose)
    std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
} else {
  success = true;
  if (proc_verbose)
    std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
}
}
TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef EPETRA_MPI
MPI_Finalize();
#endif

return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
