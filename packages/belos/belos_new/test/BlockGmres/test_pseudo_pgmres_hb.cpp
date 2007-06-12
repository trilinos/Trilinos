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
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "createEpetraProblem.hpp"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

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
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;

  bool verbose = false, proc_verbose = false;
  bool leftprec = true;  // use left preconditioning to solve these linear systems
  int frequency = -1;    // how often residuals are printed by solver
  int init_numrhs = 5;   // how many right-hand sides get solved first
  int aug_numrhs = 10;   // how many right-hand sides are augmented to the first group
  int maxrestarts = 15;  // number of restarts allowed 
  int length = 25;
  int blocksize = 5;     // blocksize used by pseudo-block GMRES
  int maxiters = -1;     // maximum iterations allowed
  std::string filename("orsirr1.hb");
  MT tol = 1.0e-5;       // relative residual tolerance
  MT aug_tol = 1.0e-5;   // relative residual tolerance for augmented systems

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("aug-tol",&aug_tol,"Relative residual tolerance used by GMRES solver for augmented systems.");
  cmdp.setOption("init-num-rhs",&init_numrhs,"Number of right-hand sides to be initially solved for.");
  cmdp.setOption("aug-num-rhs",&aug_numrhs,"Number of right-hand sides augmenting the initial solve.");
  cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
  cmdp.setOption("block-size",&blocksize,"Block size used by GMRES.");  
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
  RefCountPtr<Epetra_CrsMatrix> A;
  int return_val =Belos::createEpetraProblem(filename,NULL,&A,NULL,NULL,&MyPID);
  const Epetra_Map &Map = A->RowMap();
  if(return_val != 0) return return_val;
  proc_verbose = verbose && (MyPID==0); /* Only print on zero processor */
  //
  // *****Construct the Preconditioner*****
  //
  if (proc_verbose) cout << endl << endl;
  if (proc_verbose) cout << "Constructing ILU preconditioner" << endl;
  int Lfill = 2;
  // if (argc > 2) Lfill = atoi(argv[2]);
  if (proc_verbose) cout << "Using Lfill = " << Lfill << endl;
  int Overlap = 2;
  // if (argc > 3) Overlap = atoi(argv[3]);
  if (proc_verbose) cout << "Using Level Overlap = " << Overlap << endl;
  double Athresh = 0.0;
  // if (argc > 4) Athresh = atof(argv[4]);
  if (proc_verbose) cout << "Using Absolute Threshold Value of " << Athresh << endl;
  double Rthresh = 1.0;
  // if (argc >5) Rthresh = atof(argv[5]);
  if (proc_verbose) cout << "Using Relative Threshold Value of " << Rthresh << endl;
  //
  Teuchos::RefCountPtr<Ifpack_IlukGraph> ilukGraph;
  Teuchos::RefCountPtr<Ifpack_CrsRiluk> ilukFactors;
  //
  if (Lfill > -1) {
    ilukGraph = Teuchos::rcp(new Ifpack_IlukGraph(A->Graph(), Lfill, Overlap));
    assert(ilukGraph->ConstructFilledGraph()==0);
    ilukFactors = Teuchos::rcp(new Ifpack_CrsRiluk(*ilukGraph));
    int initerr = ilukFactors->InitValues(*A);
    if (initerr != 0) cout << "InitValues error = " << initerr;
    assert(ilukFactors->Factor() == 0);
  }
  //
  bool transA = false;
  double Cond_Est;
  ilukFactors->Condest(transA, Cond_Est);
  if (proc_verbose) {
    cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
    cout << endl;
  }

  //
  // Create the Belos preconditioned operator from the Ifpack preconditioner.
  // NOTE:  This is necessary because Belos expects an operator to apply the
  //        preconditioner with Apply() NOT ApplyInverse().
  RefCountPtr<Belos::EpetraPrecOp> Prec = rcp( new Belos::EpetraPrecOp( ilukFactors ) );

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = Map.NumGlobalElements();
  if (maxiters = -1)
    maxiters = NumGlobalElements - 1; // maximum number of iterations to run
  //
  ParameterList belosList;
  belosList.set( "Num Blocks", length );                 // Maximum number of blocks in Krylov factorization
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
  // *****Construct solution vector and random right-hand-sides *****
  //
  RefCountPtr<Epetra_MultiVector> X = rcp( new Epetra_MultiVector(Map, init_numrhs) );
  X->PutScalar( 0.0 );
  RefCountPtr<Epetra_MultiVector> B = rcp( new Epetra_MultiVector(Map, init_numrhs) );
  B->Random();
  Belos::LinearProblem<double,MV,OP> problem( A, X, B );
  if (leftprec)
    problem.setLeftPrec( Prec );
  else
    problem.setRightPrec( Prec );
  
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
    return -1;
  }
  //
  // *******************************************************************
  // *************Start the block Gmres iteration*************************
  // *******************************************************************
  //
  Teuchos::RefCountPtr< Belos::SolverManager<double,MV,OP> > solver
    = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>( rcp(&problem,false), belosList ) );
  //
  // Perform solve
  //
  Belos::ReturnType ret = solver->solve();
  
  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    cout << endl << endl;
    cout << "Dimension of matrix: " << NumGlobalElements << endl;
    cout << "Number of initial right-hand sides: " << init_numrhs << endl;
    cout << "Number of augmented right-hand sides: " << aug_numrhs << endl;
    cout << "Number of restarts allowed: " << maxrestarts << endl;
    cout << "Length of block Arnoldi factorization: " << length <<endl;
    cout << "Max number of Gmres iterations: " << maxiters << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    if (aug_tol != tol)
      cout << "Relative residual tolerance for augmented systems: " << aug_tol << endl;
    cout << endl;
  }
  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<double> actual_resids( init_numrhs );
  std::vector<double> rhs_norm( init_numrhs );
  Epetra_MultiVector R(Map, init_numrhs);
  OPT::Apply( *A, *X, R );
  MVT::MvAddMv( -1.0, R, 1.0, *B, R ); 
  MVT::MvNorm( R, &actual_resids );
  MVT::MvNorm( *B, &rhs_norm );
  if (proc_verbose) {
    cout<< "---------- Actual Residuals (normalized) ----------"<<endl<<endl;
    for ( int i=0; i<init_numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      cout<<"Problem "<<i<<" : \t"<< actRes <<endl;
      if (actRes > tol ) badRes = true;
    }
  }

  if (ret!=Belos::Converged || badRes==true) {
    if (proc_verbose)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (proc_verbose)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
  //
} // end test_bl_pgmres_hb.cpp
