// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

// The Trilinos package Galeri has many example problems.
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include "Ifpack.h"

// Include selected communicator class required by Epetra objects
#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // EPETRA_MPI

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

int
main (int argc, char *argv[])
{
  int MyPID = 0;

  // Belos solvers have the following template parameters:
  //
  //   - Scalar: The type of dot product results.
  //   - MV: The type of (multi)vectors.
  //   - OP: The type of operators (functions from multivector to
  //     multivector).  A matrix (like Epetra_CrsMatrix) is an example
  //     of an operator; an Ifpack preconditioner is another example.
  //
  // Here, Scalar is double, MV is Epetra_MultiVector, and OP is
  // Epetra_Operator.
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

#ifdef EPETRA_MPI
  MPI_Init (&argc, &argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif // EPETRA_MPI

bool verbose = false;
bool success = true;
try {
  bool proc_verbose = false;
  bool debug = false;
  bool leftprec = true;      // left preconditioning or right.
  bool pseudo = false;       // use pseudo-block or block gmres
  int frequency = -1;        // frequency of status test output
  int blocksize = 1;         // blocksize
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  int maxsubspace = 50;      // maximum number of blocks the solver can use for the subspace
  int maxrestarts = 15;      // number of restarts allowed
  int nx = 10;               // number of discretization points in each direction
  MT tol = 1.0e-5;           // relative residual tolerance
  std::string ortho = "DGKS";// orthogonalization method

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nondebug",&debug,"Print debugging information from solver.");
  cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
  cmdp.setOption("pseudo","block",&pseudo,"Use pseudo-block or block GMRES solver.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blocksize,"Block size used by GMRES.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
  cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of blocks the solver can use for the subspace.");
  cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
  cmdp.setOption("nx",&nx,"Number of discretization points in each direction of 3D Laplacian.");
  cmdp.setOption("ortho",&ortho,"Orthogonalization being used by GMRES solver.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose

  //
  // Set up the test problem.
  //
  // We use Trilinos' Galeri package to construct a test problem.
  // Here, we use a discretization of the 2-D Laplacian operator.
  // The global mesh size is nx * nx.
  //
  Teuchos::ParameterList GaleriList;
  GaleriList.set ("n", nx * nx * nx);
  GaleriList.set ("nx", nx);
  GaleriList.set ("ny", nx);
  GaleriList.set ("nz", nx);
  RCP<Epetra_Map> Map = rcp (Galeri::CreateMap ("Linear", Comm, GaleriList));
  RCP<Epetra_RowMatrix> A =
    rcp (Galeri::CreateCrsMatrix ("Laplace3D", &*Map, GaleriList));

  proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

  // Create RHS using random solution vector
  RCP<MV> B = rcp (new MV (*Map, numrhs));
  RCP<MV> X = rcp (new MV (*Map, numrhs));
  RCP<MV> Xexact = rcp (new MV (*Map, numrhs));
  Xexact->Random ();

  A->Apply( *Xexact, *B );

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
  ifpackList.set("fact: level-of-fill", 1);
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
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  //
  ParameterList belosList;
  belosList.set( "Num Blocks", maxsubspace);             // Maximum number of blocks in Krylov factorization
  belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  belosList.set( "Orthogonalization", ortho );           // Orthogonalization used by iterative solver
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
  if (leftprec) {
    problem.setLeftPrec( belosPrec );
  }
  else {
    problem.setRightPrec( belosPrec );
  }

  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  //
  //
  // *******************************************************************
  // *************Start the block Gmres iteration*************************
  // *******************************************************************
  //
  // Create an iterative solver manager.
  RCP< Belos::SolverManager<double,MV,OP> > newSolver;
  if (pseudo && (blocksize == 1))
    newSolver = rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));
  else
    newSolver = rcp( new Belos::BlockGmresSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Max number of restarts allowed: " << maxrestarts << std::endl;
    std::cout << "Max number of Gmres iterations per linear system: " << maxiters << std::endl;
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
    std::cout << "End Result: TEST FAILED" << std::endl;
} else {
  if (proc_verbose)
    std::cout << "End Result: TEST PASSED" << std::endl;
}
}
TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef EPETRA_MPI
MPI_Finalize();
#endif

return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
