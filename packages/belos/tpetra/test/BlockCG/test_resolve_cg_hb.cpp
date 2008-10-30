// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side from the problem is being used instead of multiple
// random right-hand-sides.  The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "createEpetraProblem.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
  //
#ifdef EPETRA_MPI	
  MPI_Init(&argc,&argv);
  Belos::MPIFinalize mpiFinalize; // Will call finalize with *any* return
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

  bool verbose = false, proc_verbose = false;
  bool pseudo = false;   // use pseudo block GMRES to solve this linear system.
  int frequency = -1;
  int blocksize = 1;
  int numrhs = 1;
  int maxrestarts = 15; // number of restarts allowed 
  int maxiters = -1;    // maximum number of iterations allowed per linear system
  std::string filename("bcsstk14.hb");
  MT tol = 1.0e-5;  // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("pseudo","regular",&pseudo,"Use pseudo-block GMRES to solve the linear systems.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
  cmdp.setOption("blocksize",&blocksize,"Block size used by GMRES.");
  cmdp.setOption("maxiters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
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
  int return_val =Belos::createEpetraProblem(filename,NULL,&A,&B,&X,&MyPID);
  if(return_val != 0) return return_val;
  const Epetra_Map &Map = A->RowMap();
  proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  //
  ParameterList belosList;
  belosList.set( "Num Blocks", maxiters );               // Maximum number of blocks in Krylov factorization
  belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
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
  // *************Start the block CG iteration*************************
  // *******************************************************************
  //
  Teuchos::RCP< Belos::SolverManager<double,MV,OP> > solver;
  solver = Teuchos::rcp( new Belos::BlockCGSolMgr<double,MV,OP>( rcp(&problem,false), rcp(&belosList,false) ) );
  //
  // Perform solve
  //
  Belos::ReturnType ret = solver->solve();
  //
  // Get the number of iterations for this solve.
  //
  int numIters = solver->getNumIters();
  std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;

  if (ret!=Belos::Converged) {
    if (proc_verbose)
      std::cout << "End Result: TEST FAILED" << std::endl;	
    return -1;
  }
  //
  //
  // Compute actual residuals.
  //
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(Map, numrhs);
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
  //
  //
  // Construct a second unpreconditioned linear problem instance.
  //
  RCP<Epetra_MultiVector> X2 = MVT::Clone(*X, numrhs); 
  MVT::MvInit( *X2, 0.0 );
  Belos::LinearProblem<double,MV,OP> problem2( A, X2, B );
  problem2.setLabel("Belos Resolve");
  set = problem2.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  //
  // *******************************************************************
  // *************Start the block CG iteration*************************
  // *******************************************************************
  //
  // Create the solver without either the problem or parameter list.
  solver = Teuchos::rcp( new Belos::BlockCGSolMgr<double,MV,OP>() );

  // Get the valid list of parameters from the solver and print it.
  RCP<const Teuchos::ParameterList> validList = solver->getValidParameters();
  if (proc_verbose) {
    std::cout << std::endl << "Valid parameters from the block CG solver manager:" << std::endl;
    std::cout << *validList << std::endl;
  }
  //
  // Set the problem after the solver construction.
  solver->setProblem( rcp( &problem2, false ) );

  // Set the parameter list after the solver construction.
  belosList.set( "Timer Label", "Belos Resolve" );         // Set timer label to discern between the two solvers.
  solver->setParameters( rcp( &belosList, false ) );
  //
  // Perform solve
  //
  ret = solver->solve();
  //
  // Get the number of iterations for this solve.
  //
  numIters = solver->getNumIters();
  std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;

  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<double> actual_resids2( numrhs );
  Epetra_MultiVector resid2(Map, numrhs);
  OPT::Apply( *A, *X2, resid2 );
  MVT::MvAddMv( -1.0, resid2, 1.0, *B, resid2 ); 
  MVT::MvNorm( resid2, actual_resids2 );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Residuals 2 (normalized) ----------"<<std::endl<<std::endl;
    for ( int i=0; i<numrhs; i++) {
      std::cout<<"Problem "<<i<<" : \t"<< actual_resids2[i]/rhs_norm[i] <<std::endl;
      if ( ( actual_resids2[i] - actual_resids[i] ) > SCT::prec() ) { 
        badRes = true;
      }
    }
  }
  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }

  if (ret!=Belos::Converged || badRes) {
    if (proc_verbose)
      std::cout << "End Result: TEST FAILED" << std::endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (proc_verbose)
    std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
  //
} // end test_resolve_gmres_hb.cpp
