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
// This driver reads a problem from a file, which can be in Harwell-Boeing (*.hb),
// Matrix Market (*.mtx), or triplet format (*.triU, *.triS).  The right-hand side
// from the problem, if it exists, will be used instead of multiple random
// right-hand-sides.  The initial guesses are all set to zero.  An ILU preconditioner
// is constructed using the Ifpack factory.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosGCRODRSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

#include "Ifpack.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

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

  bool verbose = false, debug = false, proc_verbose = false;
  bool leftprec = true;      // left preconditioning or right.
  int frequency = -1;        // frequency of status test output.
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  int maxsubspace = 250;     // maximum number of blocks the solver can use for the subspace
  int recycle = 50;          // maximum size of recycle space
  int maxrestarts = 15;      // maximum number of restarts allowed 
  std::string filename("orsirr1.hb");
  std::string ortho("IMGS");
  MT tol = 1.0e-10;          // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nondebug",&debug, "Print debugging information from solver.");
  cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
  cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of blocks the solver can use for the subspace.");
  cmdp.setOption("recycle",&recycle,"Number of vectors in recycle space.");
  cmdp.setOption("max-cycles",&maxrestarts,"Maximum number of cycles allowed for GCRO-DR solver.");
  cmdp.setOption("ortho-type",&ortho,"Orthogonalization type. Must be one of DGKS, ICGS, IMGS.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose

  //
  // *************Get the problem*********************
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
    X->Seed();
    X->Random();
    OPT::Apply( *A, *X, *B );
    X->PutScalar( 0.0 );
  }
  else {
    X = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecX);
    B = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecB);
  }
  //
  // ************Construct preconditioner*************
  //
  ParameterList ifpackList;

  // allocates an IFPACK factory. No data is associated
  // to this object (only method Create()).
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  std::string PrecType = "ILU"; // incomplete LU
  int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.

  RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create(PrecType, &*A, OverlapLevel) );
  assert(Prec != Teuchos::null);

  // specify parameters for ILU
  ifpackList.set("fact: drop tolerance", 1e-9);
  ifpackList.set("fact: ilut level-of-fill", 1.0);
  // the combine mode is on the following:
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  // Their meaning is as defined in file Epetra_CombineMode.h
  ifpackList.set("schwarz: combine mode", "Add");
  // sets the parameters
  IFPACK_CHK_ERR(Prec->SetParameters(ifpackList));

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Builds the preconditioners, by looking for the values of
  // the matrix.
  IFPACK_CHK_ERR(Prec->Compute());

  // Create the Belos preconditioned operator from the Ifpack preconditioner.
  // NOTE:  This is necessary because Belos expects an operator to apply the
  //        preconditioner with Apply() NOT ApplyInverse().
  RCP<Belos::EpetraPrecOp> belosPrec = rcp( new Belos::EpetraPrecOp( Prec ) );

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements - 1; // maximum number of iterations to run
  //
  ParameterList belosList;
  belosList.set( "Num Blocks", maxsubspace );            // Maximum number of blocks in Krylov factorization
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  belosList.set( "Num Recycled Blocks", recycle );       // Number of vectors in recycle space
  belosList.set( "Orthogonalization", ortho );           // Orthogonalization type
  if (numrhs > 1) {
    belosList.set( "Show Maximum Residual Norm Only", true );  // Show only the maximum residual norm
  }
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
  // *******Construct a preconditioned linear problem********
  //
  RCP<Belos::LinearProblem<double,MV,OP> > problem
    = rcp( new Belos::LinearProblem<double,MV,OP>( A, X, B ) );
  if (leftprec) {
    problem->setLeftPrec( belosPrec );
  }
  else {
    problem->setRightPrec( belosPrec );
  }    
  bool set = problem->setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  
  // Create an iterative solver manager.
  RCP< Belos::SolverManager<double,MV,OP> > solver
    = rcp( new Belos::GCRODRSolMgr<double,MV,OP>(problem, rcp(&belosList,false)));
  
  //
  // *******************************************************************
  // *************Start the block Gmres iteration*************************
  // *******************************************************************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Number of restarts allowed: " << maxrestarts << std::endl;
    std::cout << "Max number of Gmres iterations per restart cycle: " << maxiters << std::endl; 
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }
  //
  // Perform solve
  //
  Belos::ReturnType ret = solver->solve();
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

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

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
}
